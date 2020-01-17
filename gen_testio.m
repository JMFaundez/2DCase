%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Generate test for I/O (.u)                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/scratch/nicolo/work/matlab-library/nek/')
clc, close all, clear all



%% Parameters

% Grid
gridname = 'wing';
nelx = 400;
nely = 50;

% Where?
%    x    y    z
P = [ .50 -.05 0.  ];

% What?
%    u    v    w
U = [1.   0.   0.  ;
     0.   1.   0.  ];

% Output files
nekfile{1} = 'test-u.u';
nekfile{2} = 'test-v.u';
hptsfile = 'hpts.in';



%% Get GLL points

% read flow field
[nekdata,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = ...
                          readnek([gridname,'-init/',gridname,'0.f00001']);
N = lr1(1);

% reshape gll points
[xx,yy] = reshapenek(nekdata,nelx,nely);

clear nekdata



%% Closest point
dd = sqrt((xx-P(1)).^2 + (yy-P(2)).^2);
[dQ,jQ] = min(dd); [dQ,iQ] = min(dQ); jQ = jQ(iQ);

Q = [xx(jQ,iQ) yy(jQ,iQ) 0]; dQ



%% Sensor field
for i = 1:size(U,1)
    uu = zeros(size(xx));
    vv = zeros(size(xx));
    
    uu(jQ,iQ) = U(1,i);
    vv(jQ,iQ) = U(2,i);
    
    % write file
    meshdata = zeros(size(xx,1),size(xx,2),2+2+1+1);
    meshdata(:,:,1) = xx;
    meshdata(:,:,2) = yy;
    meshdata(:,:,3) = uu;
    meshdata(:,:,4) = vv;
    
    nekdata = demeshnek(meshdata,N);
    writenek(nekfile{i},nekdata,lr1,elmap,time,istep,fields,emode,wdsz,etag);
 
    clear nekdata meshdata
end



%% hpts file
fid = fopen(hptsfile,'w');
fprintf(fid,' %d\n',1);
fprintf(fid,' %23.16E %23.16E %23.16E\n',Q(1),Q(2),Q(3));
fid = fclose(fid);



%% Plot
figure(1); clf; hold on
xxel = xx(1:N-1:end,1:N-1:end);
yyel = yy(1:N-1:end,1:N-1:end);
zzel = xx(1:N-1:end,1:N-1:end) * 0;
mesh(xxel,yyel,zzel,'EdgeColor','k');
plot(P(1),P(2),'or',Q(1),Q(2),'+b','LineWidth',2,'MarkerSize',10);
legend('element mesh','P','Q')
axis equal; grid on
xlabel('x/c'); ylabel('y/c');