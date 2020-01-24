function gen_rea_base()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Generate NEK input files (.rea)                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all, clear all

write_file = 1;
plot_stuff = 0;
curved_el = 0;

% -- number of points for the grid that will be mapped
nelx = 200; % along the profile
nely =  40; % normal to the profile
n = 700;
mesh_number = 2;

% Grid
gridname = ['mesh_',num2str(mesh_number)];

% uncomment the following line if nelx or nely changes
gen_gri(gridname, nelx, nely, n, mesh_number);


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
if plot_stuff
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
end

%% Build mesh

% gridpoints numbering
iind = 1:(nelx+1)*(nely+1);
iind = reshape(iind,nely+1,nelx+1);

% generate elements
nel = nelx*nely;
iel = 0;
iel_c = 1;

% Load the midpoints of the curved elements
load(['mid_points_',num2str(mesh_number),'.mat'], 'midpoints')
EL = struct('nodenum',[],'nodes',[]);
Ec = struct('edge',[],'El',[], 'C',[], 'type',[]);
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
        % The only BC wall is at the airfoil profile and the first edge
        % represents the curved face
        if EL(iel).BC(1) == 'W' && curved_el
          Ec(iel_c).edge = 1;
          Ec(iel_c).El = iel;
          Ec(iel_c).C = [midpoints.x(iel_c) midpoints.y(iel_c) 0 0 0];
          Ec(iel_c).type = 'm';
          iel_c = iel_c + 1;
        end
    end
end
%Store the elements 
save(['EL_',num2str(mesh_number),'.mat'],'EL') 
save(['Ec_',num2str(mesh_number),'.mat'], 'Ec')

if write_file
    system(['rm ','GLL/',gridname,'-rea.rea']);
    write_rea(EL,Ec,['GLL/',gridname,'-rea'],2);
    system(['rm ','base/',gridname,'-rea.rea']);
    write_rea(EL,Ec,['base/',gridname,'-rea'],2);
    system(['rm ','fringe/',gridname,'-rea.rea']);
    write_rea(EL,Ec,['fringe/',gridname,'-rea'],2);
end
clc, clear all
