function fid = write_grid(name,xx,yy,ii,uu,vv,pp,fr)
%   
%   fid = write_grid(name,xx,yy,ii,uu,vv,pp,rr)
%   

fid = fopen([name,'.grid'],'w+');

% reshape data
[ny,nx] = size(xx);

x = reshape(xx,1,nx*ny);
y = reshape(yy,1,nx*ny);
i = reshape(ii,1,nx*ny,4);
u = reshape(uu,1,nx*ny);
v = reshape(vv,1,nx*ny);
p = reshape(pp,1,nx*ny);
f = reshape(fr,1,nx*ny);

% header
fprintf(fid,'# nx       ny      \n %8i %8i\n',nx,ny);

% data
fprintf(fid,['# x                     ',...
              ' y                     ',...
              ' boundary-id N'         ,...
              ' boundary-id E'         ,...
              ' boundary-id S'         ,...
              ' boundary-id W'         ,...
              ' u                     ',...
              ' v                     ',...
              ' p                     ',...
              ' fringe                ','\n']);
fprintf(fid,' %22.15e %22.15e %10i      %12i      %10i      %10i      %22.15e %22.15e %22.15e %22.15e\n',...
            [ x     ; y     ; i(:,:,1); i(:,:,2); i(:,:,3); i(:,:,4); u     ; v     ; p     ; f     ]);

fid = fclose(fid);