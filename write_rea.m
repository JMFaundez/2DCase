function fid = write_rea(EL, Ec,name,dim, varargin)
% This function reads Elements and curved elements to create the rea file
% This rea file does not include the parameters since it is just used 
%  to create a .re2 file by using the tool reatore2 from Nek.
%
% EL : Elements with their nodes and nodes coordinates
% Ec : Curved elements following the Nek format
% name: name of the rea file that will be stores
% dim: number of dimensions

bcflag = 0;

% -read varargin
for i = 2:2:length(varargin)
    switch varargin{i-1}
        case 'BCfile'
            bcflag = 1;
            bcname = varargin{i};
    end
end

% Dimension check
if dim ~= 2
    disp('Error: only 2D implemented');
    fid = -1; return
end

nel = length(EL);
nec = length(Ec);

% Open file
fid = fopen([name,'.rea'],'w+');
% Grid
fprintf(fid,' **MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.\n');
fprintf(fid,' %11i %2i %11i           NEL,NDIM,NELV\n',nel,2,nel);
for iel = 1:nel
    fprintf(fid,'            ELEMENT %11i [    1a]  GROUP  0\n',iel);
    fprintf(fid,'  %22.15E %22.15E %22.15E %22.15E\n',EL(iel).nodes(1,:));
    fprintf(fid,'  %22.15E %22.15E %22.15E %22.15E\n',EL(iel).nodes(2,:));
end


% Curved edges
fprintf(fid,'  ***** CURVED SIDE DATA *****\n');
fprintf(fid,'%5i Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n', nec);
for iec = 1:nec
    if nel<1000
        fprintf(fid,'%3i%3i%14.6E%14.6E%14.6E%14.6E%14.6E %s\n', Ec(iec).edge,Ec(iec).El,Ec(iec).C(1),Ec(iec).C(2),Ec(iec).C(3),Ec(iec).C(4),Ec(iec).C(5), Ec(iec).type);
    elseif nel< 1000000
        fprintf(fid,'%2i%6i%14.6E%14.6E%14.6E%14.6E%14.6E %s\n', Ec(iec).edge,Ec(iec).El,Ec(iec).C(1),Ec(iec).C(2),Ec(iec).C(3),Ec(iec).C(4),Ec(iec).C(5), Ec(iec).type);
    else
        fprintf(fid,'%2i%12i%14.6E%14.6E%14.6E%14.6E%14.6E %s\n', Ec(iec).edge,Ec(iec).El,Ec(iec).C(1),Ec(iec).C(2),Ec(iec).C(3),Ec(iec).C(4),Ec(iec).C(5), Ec(iec).type);
    end
end

% Boundary conditions
fprintf(fid,'  ***** BOUNDARY CONDITIONS *****\n');
fprintf(fid,'  ***** FLUID BOUNDARY CONDITIONS *****\n');

if nel < 1000000
    for iel = 1:nel
        for ind = 1:4
            fprintf(fid,' %s   %5i%1i   0.00000       0.00000       0.00000       0.00000       0.00000\n',...
                    EL(iel).BC(ind),iel,ind);
        end
    end
end

fprintf(fid,'  ***** NO THERMAL BOUNDARY CONDITIONS *****\n');
if bcflag
    fprintf(fid,'   1 PRESOLVE/RESTART OPTIONS  *****\n');
    fprintf(fid,' %s XUPT\n',bcname);
else
    fprintf(fid,'   0 PRESOLVE/RESTART OPTIONS  *****\n');
end
    
% other parameters
fprintf(fid,'   7         INITIAL CONDITIONS *****\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'C Default\n');
fprintf(fid,'  ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q\n');
fprintf(fid,'   4                 Lines of Drive force data follow\n');
fprintf(fid,'C\n');
fprintf(fid,'C\n');
fprintf(fid,'C\n');
fprintf(fid,'C\n');
fprintf(fid,'  ***** Variable Property Data ***** Overrrides Parameter data.\n');
fprintf(fid,'   1 Lines follow.\n');
fprintf(fid,'   0 PACKETS OF DATA FOLLOW\n');
fprintf(fid,'  ***** HISTORY AND INTEGRAL DATA *****\n');
fprintf(fid,'   0   POINTS.  Hcode, I,J,H,IEL\n');
fprintf(fid,'  ***** OUTPUT FIELD SPECIFICATION *****\n');
fprintf(fid,'   6 SPECIFICATIONS FOLLOW\n');
fprintf(fid,'   T      COORDINATES\n');
fprintf(fid,'   T      VELOCITY\n');
fprintf(fid,'   T      PRESSURE\n');
fprintf(fid,'   T      TEMPERATURE\n');
fprintf(fid,'   F      TEMPERATURE GRADIENT\n');
fprintf(fid,'   0      PASSIVE SCALARS\n');
fprintf(fid,'  ***** OBJECT SPECIFICATION *****\n');
fprintf(fid,'       0 Surface Objects\n');
fprintf(fid,'       0 Volume  Objects\n');
fprintf(fid,'       0 Edge    Objects\n');
fprintf(fid,'       0 Point   Objects\n');
      
% close file
fid = fclose(fid);
