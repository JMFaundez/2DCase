%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Generate NEK input files (.rea, .bc)                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/scratch/josfa/matlab-tools/nek/')
clc, close all, clear all

write_file = 1;

% Grid
gridname = 'FST_naca0008'; runnek = 0; % -> run Nek in gridname-init by yourself because you need the right values in the SIZE file
% -- number of points for the grid that will be mapped
nelx = 50; % along the profile
nely =  10; % normal to the profile
n = 600;
gen_gri(gridname, nelx, nely, n);

% -- boundary conditions for each boundary (W: wall, v: Dirichlet, O: Neumann)
bc{0+1} = 'E'; % no bc for internal nodes
bc{1+1} = 'W';
bc{2+1} = 'o';
bc{3+1} = 'v';
bc{4+1} = 'v';
bc{5+1} = 'v';
bc{6+1} = 'o';


%% Load grid
[xx,yy,ii,uu,vv,pp,fr] = read_grid(gridname); [nely,nelx] = size(xx(2:end,2:end));

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
save('EL.mat','EL')

if write_file
system(['rm ',gridname,'-init/',gridname,'-rea.rea']);
% system(['mkdir ',gridname,'-init']);
write_rea(EL,[gridname,'-init/',gridname,'-rea'],2);
end

