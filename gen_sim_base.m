%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Generate NEK input files (.rea, .bc)                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath('/scratch/nicolo/work/matlab-library/nek/')
addpath('/scratch/josfa/matlab-tools/nek/')
clc, close all, clear all



%% Parameters

% Grid
gridname = 'FST_naca0008'; runnek = 0; % -> run Nek in gridname-init by yourself because you need the right values in the SIZE file


% Initial condition (if none data from grid file)
% icdata = 'bc_morino/alphaexp-2.5.mat';
% bldata = 'bc_morino/bl-data.mat';


% Simulation
simname = 'FST_naca0008-base';

Re = 5.333333e5;%3.75e6;
%U  = 10; <- needed only if icdata exists

% -- boundary conditions for each boundary (W: wall, v: Dirichlet, O: Neumann)
bc{0+1} = 'E'; % no bc for internal nodes
bc{1+1} = 'W';
bc{2+1} = 'o';
bc{3+1} = 'v';
bc{4+1} = 'v';
bc{5+1} = 'v';
bc{6+1} = 'o';

% -- fringe
stfr = 0.0;



%% Load grid
[xx,yy,ii,uu,vv,pp,fr] = read_grid(gridname); [nely,nelx] = size(xx(2:end,2:end));



%% Initial/Boundary condition 
if exist('icdata','var')
    
    % generate from potential solution
    mo = load(icdata);
    [uu,vv,pp] = morinoflowfield(mo.pan+mo.dpan,mo.cll+mo.dcll,mo.alpha,mo.phi,xx,yy);
    
    % rescale velocity
    uu = U * uu; vv = U * vv; pp = 1 - (uu.^2 + vv.^2)/2;
    
    % fix boundary layer
    bl = load(bldata);
    dstfun = @(m) interp1(bl.m,bl.dst,m,'spline');
    ublfun = @(eta,m) interp2(bl.m,bl.yy,squeeze(bl.ff(2,:,:)),m,eta,'spline',1);
    etamax = min([max(bl.yy) 20])
    
    % wall normal direction
    n = [xx(end,:) - xx(1,:);
         yy(end,:) - yy(1,:)]; n = n./([1 1]' * sqrt(n(1,:).^2 + n(2,:).^2));
    
    % tangential coordinate
    dx = xx(1,2:end) - xx(1,1:end-1);
    dy = yy(1,2:end) - yy(1,1:end-1);

    ds = sqrt(dx.^2+dy.^2);
    spr = cumsum([0 ds]); [~,ispr0] = min(abs(xx(1,:))); spr = spr - spr(ispr0);
    
    % boundary layer characteristics
    mbl = interp1(mo.scll,mo.m,spr,'spline');
    deltabl = interp1(mo.scll,mo.delta,spr,'spline');
    
    for i = 1:nelx+1
        
        % boundary layer coordinate
        etabl = n(:,i)' * [ xx(:,i)' - xx(1,i) ;
                            yy(:,i)' - yy(1,i) ] ./ deltabl(i);
        jj = etabl <= etamax;
        
        % regularize velocity (close to the wall)
        method = 'linear';
        uu(jj,i) = interp1(etabl(~jj),uu(~jj,i),etabl(jj),method,'extrap');
        vv(jj,i) = interp1(etabl(~jj),vv(~jj,i),etabl(jj),method,'extrap');
        pp(jj,i) = 1 - (uu(jj,i).^2 + vv(jj,i).^2)/2;
        
        % map boundary layer to grid
        uu(jj,i) = ublfun(etabl(jj),mbl(i)) .* uu(jj,i);
        vv(jj,i) = ublfun(etabl(jj),mbl(i)) .* vv(jj,i);
        
    end
    
end



%% Plot grid and ic/bc
figure(1); clf; hold on; set (1,'Units','normalized','Position',[.5 .5 .5 .5]);

%hh = surf(xx,yy,uu,'EdgeColor','none'); hc = colorbar('NO'); xlabel(hc,'u')
%hh = surf(xx,yy,vv,'EdgeColor','none'); hc = colorbar('NO'); xlabel(hc,'v')
hh = surf(xx,yy,pp,'EdgeColor','none'); hc = colorbar('NO'); xlabel(hc,'p')
os = max(max(get(hh,'ZData')));

clrs = 'rgbcmy';
for i = 1:max(max(max(ii)))
    bd = sum(ii==i,3)~=0;
    hbd(i) = plot3(xx(bd),yy(bd),os+0*xx(bd),...
                   clrs(mod(i-1,6)+1),'LineWidth',1.5); % boundaries
    lgd{i} = ['boundary ',(num2str(i)),' (',bc{i+1},')'];
end
        
legend(hbd,lgd,'Location','SO','Orientation','Horizontal');

hold off; view(2); axis image; grid on
xlabel('x/c'); ylabel('y/c'); title('Grid'); drawnow



%% Build mesh

% gridpoints numbering
iind = 1:(nelx+1)*(nely+1);
iind = reshape(iind,nely+1,nelx+1);

% generate elements
nel = nelx*nely;

iel = 0;
for i = 1:nelx
    for j = 1:nely
        % element id
        iel = iel+1;
        
        % corners
        EL(iel).nodenum = [iel iind(j,i) iind(j,i+1) iind(j+1,i+1) iind(j+1,i)];
        EL(iel).nodes(1,:) = [xx(j,i) xx(j,i+1) xx(j+1,i+1) xx(j+1,i)];
        EL(iel).nodes(2,:) = [yy(j,i) yy(j,i+1) yy(j+1,i+1) yy(j+1,i)];
        
        % bc
        EL(iel).BC = [ bc{ii(j  ,i  ,1)+1} ...
                       bc{ii(j  ,i+1,4)+1} ...
                       bc{ii(j+1,i+1,3)+1} ...
                       bc{ii(j+1,i  ,2)+1} ];
        
    end
end

system(['rm ',gridname,'-init/',gridname,'a.rea']);
% system(['mkdir ',gridname,'-init']);
write_rea(EL,[gridname,'-init/',gridname,'a'],2);


