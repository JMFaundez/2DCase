function [xx,yy,ii,uu,vv,pp,fr] = read_grid(name)
%   
%   [xx,yy,ii,uu,vv,pp,rr] = read_grid(name)
%   

fid = fopen([name,'.grid'],'r');

% grid size
data = textscan(fid,'',2,'HeaderLines',1); nx = data{1}; ny = data{2}; clear data

% read grid data
data = textscan(fid,'','HeaderLines',1);

x = data{1};
y = data{2};
i(:,1,1) = data{3}; i(:,1,2) = data{4}; i(:,1,3) = data{5}; i(:,1,4) = data{6};
u = data{7};
v = data{8}; 
p = data{9};
f = data{10}; 
clear data


% reshape data
xx = reshape(x,ny,nx);
yy = reshape(y,ny,nx);
ii = reshape(i,ny,nx,4);
uu = reshape(u,ny,nx);
vv = reshape(v,ny,nx);
pp = reshape(p,ny,nx);
fr = reshape(f,ny,nx);

fid = fclose(fid);