clc, close all, clear all
addpath('/scratch/josfa/matlab-tools/nek/')

mesh_n = 2;
gridname = ['mesh_',num2str(mesh_n)];

%% Load grid
[xx,yy,ii,uu,vv,pp,fr] = read_grid(gridname); [nely,nelx] = size(xx(2:end,2:end));
% -- fringe
stfr = 0 %2/(1.2e-2);
simname = ['mesh_',num2str(mesh_n),'_BC'];
% -- boundary conditions for each boundary (W: wall, v: Dirichlet, O: Neumann)
bc{0+1} = 'E'; % no bc for internal nodes
bc{1+1} = 'W';
bc{2+1} = 'o';
bc{3+1} = 'v';
bc{4+1} = 'v';
bc{5+1} = 'v';
bc{6+1} = 'o';

Re = 5.333333e5;%3.75e6;
%% Get GLL points

% read flow field
[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(['GLL/',gridname,'0.f00001']);
nel = length(elmap);
load(['EL_',num2str(mesh_n),'.mat'],'EL');
% save GLL points in EL structure
for iel = 1:nel
    EL(iel).GLL(:,1) = squeeze(nekdata(iel,:,1));
    EL(iel).GLL(:,2) = squeeze(nekdata(iel,:,2));
end

% reshape gll points
[xxic,yyic] = reshapenek(nekdata,nelx,nely);

clear nekdata




%% Interpolate data from .grid to GLL points

% create an index-based grid
% - grid
igr = 1:(lr1(1)-1):size(xxic,2);
jgr = 1:(lr1(2)-1):size(xxic,1);

[iigr,jjgr] = meshgrid(igr,jgr);

% - gll
ptmap = (1-cos(linspace(0,pi,lr1(1))))/2 * (lr1(1)-1);
igll = zeros(1,size(xxic,2));
for ielx = 1:nelx
    ii = (0:lr1(1)-1) + (lr1(1)-1)*(ielx-1) + 1;
    igll(ii) = ptmap + (lr1(1)-1)*(ielx-1) + 1;
end

ptmap = (1-cos(linspace(0,pi,lr1(2))))/2 * (lr1(2)-1);
jgll = zeros(1,size(xxic,1));
for iely = 1:nely
    jj = (0:lr1(2)-1) + (lr1(2)-1)*(iely-1) + 1;
    jgll(jj) = ptmap + (lr1(2)-1)*(iely-1) + 1;
end
    
[iigll,jjgll] = meshgrid(igll,jgll);

% interpolate fringe
method = 'cubic';
uuic = interp2(iigr,jjgr,uu,iigll,jjgll,method);
vvic = interp2(iigr,jjgr,vv,iigll,jjgll,method);
ppic = interp2(iigr,jjgr,pp,iigll,jjgll,method);
ffic = interp2(iigr,jjgr,fr,iigll,jjgll,method);
ffic(ffic<2*eps) = 0; ffic = ffic * stfr/max(max(abs(ffic)));



%% Pressure for outflow bc nodes as negative fringe (check simname.usr for consistency)

% lgl derivative
[Dlgl,xlgl] = lgldif(lr1(1)); xlgl = (1-xlgl)/2; Dlgl = -Dlgl*2;

%
if strcmp(bc{2+1},'o')
    ibc  = size(xxic,2);
    %enx = xxic(:,ibc) - xxic(:,ibc-1);
    ibel = ibc-(lr1(1)-1):ibc;
    lbel = sqrt( (xxic(:,ibel(end)) - xxic(:,ibel(1))).^2 ...
               + (yyic(:,ibel(end)) - yyic(:,ibel(1))).^2 );
           
    % ambient pressure (pa = -1/Re dU/dx + p)
    dudx = uuic(:,ibel) * Dlgl(end,:)' ./ lbel;
    pa = ppic(:,ibc) - 1/Re * dudx;
    
    % pa is negative fringe
    ffic(:,ibc) = pa;
end

if strcmp(bc{6+1},'o')
    ibc  = 1;
    ibel = ibc:ibc+(lr1(1)-1);
    lbel = sqrt( (xxic(:,ibel(end)) - xxic(:,ibel(1))).^2 ...
               + (yyic(:,ibel(end)) - yyic(:,ibel(1))).^2 );
           
    % ambient pressure (pa = p-1/Re dU/dn)
    dudx = uuic(:,ibel) * Dlgl(end,:)' ./ lbel;
    pa = ppic(:,ibc) - 1/Re * dudx;
    
    % pa is negative fringe
    ffic(:,ibc) = -pa;
end




%% Reshape data in EL structure 

for ielx = 1:nelx
	for iely = 1:nely
            
        iel = iely + nely*(ielx-1);
            
        ii = (0:lr1(1)-1) + (lr1(1)-1)*(ielx-1) + 1;
        jj = (0:lr1(2)-1) + (lr1(2)-1)*(iely-1) + 1;
        
        EL(iel).FLD(:,1) = reshape(uuic(jj,ii)',lr1(1)*lr1(2),1);
        EL(iel).FLD(:,2) = reshape(vvic(jj,ii)',lr1(1)*lr1(2),1);
        EL(iel).PRS(:,1) = reshape(ppic(jj,ii)',lr1(1)*lr1(2),1);
        EL(iel).FRI(:,1) = reshape(ffic(jj,ii)',lr1(1)*lr1(2),1);
            
	end
end

% check interpolation for NaN
if (sum(sum(isnan(uuic))) ~= 0)|(sum(sum(isnan(vvic))) ~= 0)|(sum(sum(isnan(ppic))) ~= 0)
    disp(' Error: NaN values from griddata.')
    return
end

clear xgll ygll ugll vgll pgll fgll



%% Write initial/boundary condition with fringe

% reshape data for writing
nekdata = zeros(nel,lr1(1)*lr1(2),6);
for iel = 1:nel
    nekdata(iel,:,1) = EL(iel).GLL(:,1);
    nekdata(iel,:,2) = EL(iel).GLL(:,2);
    nekdata(iel,:,3) = EL(iel).FLD(:,1);
    nekdata(iel,:,4) = EL(iel).FLD(:,2);
    nekdata(iel,:,5) = EL(iel).PRS(:,1);
    nekdata(iel,:,6) = EL(iel).FRI(:,1);
end

% write initial condition file
if stfr>0
    status =writenek(['fringe/',simname,'.bc'],nekdata,lr1,elmap,0,0,fields,emode,wdsz,etag);
    status =writenek(['fringe/',simname,'.bc0.f00001'],nekdata,lr1,elmap,0,0,fields,emode,wdsz,etag);
else
    status =writenek(['base/',simname,'.bc'],nekdata,lr1,elmap,0,0,fields,emode,wdsz,etag);  
    status =writenek(['base/',simname,'.bc0.f00001'],nekdata,lr1,elmap,0,0,fields,emode,wdsz,etag); 
end
%status = writenek(['base-torun/',simname,'.bc0.f00001'],nekdata,lr1,elmap,0,0,fields,emode,4,etag);

