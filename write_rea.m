function fid = write_rea(EL,name,dim,varargin)
%
%   fid = write_rea(EL,name,dim,varargin)
%

bcflag = 0;
param = init_param;

% -read varargin
for i = 2:2:length(varargin)
    switch varargin{i-1}
        case 'BCfile'
            bcflag = 1;
            bcname = varargin{i};
    end
end

if mod(length(varargin),2) == 1
    param = varargin{end};
    if length(param) ~= 118
        disp('Error: param has to contain 118 elements')
        fid = -2; return
    end
end

% Dimension check
if dim ~= 2
    disp('Error: only 2D implemented');
    fid = -1; return
end

nel = length(EL);


% Open file
fid = fopen([name,'.rea'],'w+');


% Header
fprintf(fid,' ****** PARAMETERS *****\n');
fprintf(fid,'   2.6000      NEKTON VERSION\n');
fprintf(fid,'   %1i DIMENSIONAL RUN\n',dim);
fprintf(fid,'         118  PARAMETERS FOLLOW\n');
% Parameters
i = 0;
i = i+1; fprintf(fid,'  % 9.5E P001: DENSITY\n'                          ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P002: VISCOS\n'                           ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P003: : : BETAG\n'                        ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P004: : : GTHETA\n'                       ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P005: : : PGRADX\n'                       ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P006:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P007: RHOCP\n'                            ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P008: CONDUCT\n'                          ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P009:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P010: FINTIME\n'                          ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P011: NSTEPS\n'                           ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P012: DT\n'                               ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P013: IOCOMM\n'                           ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P014: IOTIME\n'                           ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P015: IOSTEP\n'                           ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P016: PSSOLVER: 0=default\n'              ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P017:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P018: GRID < 0 --> # cells on screen\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P019: INTYPE\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P020: NORDER\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P021: DIVERGENCE\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P022: HELMHOLTZ\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P023: NPSCAL\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P024: TOLREL\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P025: TOLABS\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P026: COURANT/NTAU\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P027: TORDER\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P028: TORDER: mesh velocity (0: p28=p27)\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P029: = magnetic visc if > 0, = -1/Rm if < 0\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P030: > 0 ==> properties set in uservp()\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P031: NPERT: #perturbation modes\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P032: #BCs in re2 file, if > 0\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P033: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P034: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P035: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P036: : : XMAGNET\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P037: : : NGRIDS\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P038: : : NORDER2\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P039: : : NORDER3\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P040:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P041: 1-->multiplicative SEMG\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P042: 0=gmres/1=pcg\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P043: 0=semg/1=schwarz\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P044: 0=E-based/1=A-based prec.\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P045: Relaxation factor for DTFS\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P046: reserved\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P047: vnu: mesh matieral prop.\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P048: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P049: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P050: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P051:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P052: IOHIS\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P053:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P054: fixed flow rate dir: |p54|=1,2,3=x,y,z\n'                              ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P055: vol.flow rate (p54>0) or Ubar (p54<0)\n'                               ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P056: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P057: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P058:\n'                                                                     ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P059: !=0 --> full Jac. eval. for each el.\n'                                ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P060: !=0 --> init. velocity to small nonzero\n'                             ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P061:\n'                                                                     ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P062: >0 --> force byte_swap for output\n'                                   ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P063:    --> force 8-byte output\n'                                          ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P064: =1 --> perturbation restart\n'                                         ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs\n'                            ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P066: output : <0=ascii, else , nEl, binary\n'                               ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P067: restart: <0=ascii, else binary\n'                                      ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P068: iastep: freq for avg_all (0=iostep)\n'                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P069: /= 0 if restart [proper_restart.f] \n'                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P070: checkpiont dump frequency (number of time steps) [proper_restart.f]\n' ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P071: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P072: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P073:\n'                                                                     ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P074: verbose Helmholtz\n'                                                   ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P075: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P076: : :\n'                                                                 ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P077: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P078: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P079: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P080: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P081: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P082: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P083:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P084: !=0 --> sets initial timestep if p12>0\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P085: dt ratio if p84 !=0, for timesteps>0\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P086: reserved\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P087: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P088: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P089: : :   coarse grid weighting (default=10.\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P090: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P091: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P092:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P093: Number of previous pressure solns saved\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P094: start projecting velocity after p94 step\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P095: start projecting pressure after p95 step\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P096: : :   which saving algorithm 1 = discard\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P097: : :   0 == > no iterative refinement\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P098:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P099: dealiasing: <0--> off/3--> old/4--> new\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P100:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P101: Number of additional modes to filter\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P102: Dump out divergence at each time step\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P103: weight of stabilizing filter (.01)\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P104: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P105: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P106:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P107: !=0 --> add to h2 array in hlmhotz eqn\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P108: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P109: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P110: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P111: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P112: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P113: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P114: : :\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P115:\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P116: !=0: x elements for fast tensor product\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P117: !=0: y elements for fast tensor product\n'                                  ,param(i));
i = i+1; fprintf(fid,'  % 9.5E P118: !=0: z elements for fast tensor product\n'                                  ,param(i));
fprintf(fid,'      4  Lines of passive scalar data follows2 CONDUCT; 2RHOCP\n');
fprintf(fid,'   1.00000       1.00000       1.00000       1.00000       1.00000\n');
fprintf(fid,'   1.00000       1.00000       1.00000       1.00000\n');
fprintf(fid,'   1.00000       1.00000       1.00000       1.00000       1.00000\n');
fprintf(fid,'   1.00000       1.00000       1.00000       1.00000\n');
fprintf(fid,'          13   LOGICAL SWITCHES FOLLOW\n');
fprintf(fid,' T      IFFLOW\n');
fprintf(fid,' F      IFHEAT\n');
fprintf(fid,' T      IFTRAN\n');
fprintf(fid,' T F F F F F F F F F F  IFNAV & IFADVC (convection in P.S. fields)\n');
fprintf(fid,' F F T T T T T T T T T T  IFTMSH (IF mesh for this field is T mesh)\n');
fprintf(fid,' F      IFAXIS\n');
fprintf(fid,' F      IFSTRS\n');
fprintf(fid,' F      IFSPLIT\n');
fprintf(fid,' F      IFMGRID\n');
fprintf(fid,' F      IFMODEL\n');
fprintf(fid,' F      IFKEPS\n');
fprintf(fid,' F      IFMVBD\n');
fprintf(fid,' F      IFCHAR\n');
fprintf(fid,'   2.00000       2.00000      -1.00000      -1.00000     XFAC,YFAC,XZERO,YZERO\n');


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
fprintf(fid,'           0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n');


% Boundary conditions
fprintf(fid,'  ***** BOUNDARY CONDITIONS *****\n');
fprintf(fid,'  ***** FLUID BOUNDARY CONDITIONS *****\n');

if nel < 1e3
    for iel = 1:nel
        for ind = 1:4
            fprintf(fid,' %s  %3i %3i 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00\n',...
                    EL(iel).BC(ind),iel,ind);
        end
    end
elseif nel < 1e4
    for iel = 1:nel
        for ind = 1:4
            fprintf(fid,' %s   %4i%1i   0.00000       0.00000       0.00000       0.00000       0.00000\n',...
                    EL(iel).BC(ind),iel,ind);
        end
    end
else
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
