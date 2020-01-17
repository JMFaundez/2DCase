function param = init_param()
%
%   param = init_param()
%

param = zeros(118,1);

i = 0;
i = i+1; param(i) = 1.00000    ; % P001: DENSITY
i = i+1; param(i) = 1.00000    ; % P002: VISCOS
i = i+1; param(i) = 0.00000    ; % P003: : : BETAG
i = i+1; param(i) = 0.00000    ; % P004: : : GTHETA
i = i+1; param(i) = 0.00000    ; % P005: : : PGRADX
i = i+1; param(i) = 0.00000    ; % P006:
i = i+1; param(i) = 1.00000    ; % P007: RHOCP
i = i+1; param(i) = 1.00000    ; % P008: CONDUCT
i = i+1; param(i) = 0.00000    ; % P009:
i = i+1; param(i) = 0.00000    ; % P010: FINTIME
i = i+1; param(i) = 10         ; % P011: NSTEPS
i = i+1; param(i) = 5.00000E-03; % P012: DT
i = i+1; param(i) = 0.00000    ; % P013: IOCOMM
i = i+1; param(i) = 0.00000    ; % P014: IOTIME
i = i+1; param(i) = 1          ; % P015: IOSTEP
i = i+1; param(i) = 0.00000    ; % P016: PSSOLVER: 0=default
i = i+1; param(i) = 1.00000    ; % P017:
i = i+1; param(i) = 0.50000E-01; % P018: GRID < 0 --> # cells on screen
i = i+1; param(i) =-1.00000    ; % P019: INTYPE
i = i+1; param(i) = 12.0000    ; % P020: NORDER
i = i+1; param(i) = 1.00000E-09; % P021: DIVERGENCE
i = i+1; param(i) = 1.00000E-09; % P022: HELMHOLTZ
i = i+1; param(i) = 0.00000    ; % P023: NPSCAL
i = i+1; param(i) = 1.00000E-01; % P024: TOLREL
i = i+1; param(i) = 1.00000E-01; % P025: TOLABS
i = i+1; param(i) = 1.00000    ; % P026: COURANT/NTAU
i = i+1; param(i) = 3.00000    ; % P027: TORDER
i = i+1; param(i) = 0.00000    ; % P028: TORDER: mesh velocity (0: p28=p27)
i = i+1; param(i) = 0.00000    ; % P029: = magnetic visc if > 0, = -1/Rm if < 0
i = i+1; param(i) = 0.00000    ; % P030: > 0 ==> properties set in uservp()
i = i+1; param(i) = 0.00000    ; % P031: NPERT: #perturbation modes
i = i+1; param(i) = 0.00000    ; % P032: #BCs in re2 file, if > 0
i = i+1; param(i) = 0.00000    ; % P033: : :
i = i+1; param(i) = 0.00000    ; % P034: : :
i = i+1; param(i) = 0.00000    ; % P035: : :
i = i+1; param(i) = 0.00000    ; % P036: : : XMAGNET
i = i+1; param(i) = 0.00000    ; % P037: : : NGRIDS
i = i+1; param(i) = 0.00000    ; % P038: : : NORDER2
i = i+1; param(i) = 0.00000    ; % P039: : : NORDER3
i = i+1; param(i) = 0.00000    ; % P040:
i = i+1; param(i) = 0.00000    ; % P041: 1-->multiplicative SEMG
i = i+1; param(i) = 0.00000    ; % P042: 0=gmres/1=pcg
i = i+1; param(i) = 0.00000    ; % P043: 0=semg/1=schwarz
i = i+1; param(i) = 0.00000    ; % P044: 0=E-based/1=A-based prec.
i = i+1; param(i) = 0.00000    ; % P045: Relaxation factor for DTFS
i = i+1; param(i) = 0.00000    ; % P046: reserved
i = i+1; param(i) = 0.00000    ; % P047: vnu: mesh matieral prop.
i = i+1; param(i) = 0.00000    ; % P048: : :
i = i+1; param(i) = 0.00000    ; % P049: : :
i = i+1; param(i) = 0.00000    ; % P050: : :
i = i+1; param(i) = 0.00000    ; % P051:
i = i+1; param(i) = 0.00000    ; % P052: IOHIS
i = i+1; param(i) = 0.00000    ; % P053:
i = i+1; param(i) = 0.00000    ; % P054: fixed flow rate dir: |p54|=1,2,3=x,y,z
i = i+1; param(i) = 0.00000    ; % P055: vol.flow rate (p54>0) or Ubar (p54<0)
i = i+1; param(i) = 1.00000    ; % P056: : :
i = i+1; param(i) = 0.00000    ; % P057: : :
i = i+1; param(i) = 0.00000    ; % P058:
i = i+1; param(i) = 0.00000    ; % P059: !=0 --> full Jac. eval. for each el.
i = i+1; param(i) = 0.00000    ; % P060: !=0 --> init. velocity to small nonzero
i = i+1; param(i) = 0.00000    ; % P061:
i = i+1; param(i) = 0.00000    ; % P062: >0 --> force byte_swap for output
i = i+1; param(i) = 1.00000    ; % P063: =8 --> force 8-byte output\n',double);
i = i+1; param(i) = 0.00000    ; % P064: =1 --> perturbation restart
i = i+1; param(i) = 1.00000    ; % P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs
i = i+1; param(i) = 6.00000    ; % P066: output : <0=ascii, else , nEl, binary
i = i+1; param(i) = 6.00000    ; % P067: restart: <0=ascii, else binary
i = i+1; param(i) = 20000      ; % P068: iastep: freq for avg_all (0=iostep)
i = i+1; param(i) = 0.00000    ; % P069: /= 0 if restart [proper_restart.f] 
i = i+1; param(i) = 1.00000E+03; % P070: checkpiont dump frequency (number of time steps) [proper_restart.f]
i = i+1; param(i) = 0.00000    ; % P071: : :
i = i+1; param(i) = 0.00000    ; % P072: : :
i = i+1; param(i) = 0.00000    ; % P073:
i = i+1; param(i) = 0.00000    ; % P074: verbose Helmholtz
i = i+1; param(i) = 0.00000    ; % P075: : :
i = i+1; param(i) = 0.00000    ; % P076: : :
i = i+1; param(i) = 0.00000    ; % P077: : :
i = i+1; param(i) = 0.00000    ; % P078: : :
i = i+1; param(i) = 0.00000    ; % P079: : :
i = i+1; param(i) = 0.00000    ; % P080: : :
i = i+1; param(i) = 0.00000    ; % P081: : :
i = i+1; param(i) = 0.00000    ; % P082: : :
i = i+1; param(i) = 0.00000    ; % P083:
i = i+1; param(i) = 0.00000    ; % P084: !=0 --> sets initial timestep if p12>0
i = i+1; param(i) = 0.00000    ; % P085: dt ratio if p84 !=0, for timesteps>0
i = i+1; param(i) = 0.00000    ; % P086: reserved
i = i+1; param(i) = 0.00000    ; % P087: : :
i = i+1; param(i) = 0.00000    ; % P088: : :
i = i+1; param(i) = 0.00000    ; % P089: : :   coarse grid weighting (default=10.
i = i+1; param(i) = 0.00000    ; % P090: : :
i = i+1; param(i) = 0.00000    ; % P091: : :
i = i+1; param(i) = 0.00000    ; % P092:
i = i+1; param(i) = 20.0000    ; % P093: Number of previous pressure solns saved
i = i+1; param(i) = 9.00000    ; % P094: start projecting velocity after p94 step
i = i+1; param(i) = 9.00000    ; % P095: start projecting pressure after p95 step
i = i+1; param(i) = 0.00000    ; % P096: : :   which saving algorithm 1 = discard
i = i+1; param(i) = 0.00000    ; % P097: : :   0 == > no iterative refinement
i = i+1; param(i) = 0.00000    ; % P098:
i = i+1; param(i) = 3.00000    ; % P099: dealiasing: <0--> off/3--> old/4--> new
i = i+1; param(i) = 0.00000    ; % P100:
i = i+1; param(i) = 0.00000    ; % P101: Number of additional modes to filter
i = i+1; param(i) = 1.00000    ; % P102: Dump out divergence at each time step
i = i+1; param(i) = 0.01000    ; % P103: weight of stabilizing filter (.01)
i = i+1; param(i) = 0.00000    ; % P104: : :
i = i+1; param(i) = 0.00000    ; % P105: : :
i = i+1; param(i) = 0.00000    ; % P106:
i = i+1; param(i) = 0.00000    ; % P107: !=0 --> add to h2 array in hlmhotz eqn
i = i+1; param(i) = 0.00000    ; % P108: : :
i = i+1; param(i) = 0.00000    ; % P109: : :
i = i+1; param(i) = 0.00000    ; % P110: : :
i = i+1; param(i) = 0.00000    ; % P111: : :
i = i+1; param(i) = 0.00000    ; % P112: : :
i = i+1; param(i) = 0.00000    ; % P113: : :
i = i+1; param(i) = 0.00000    ; % P114: : :
i = i+1; param(i) = 0.00000    ; % P115:
i = i+1; param(i) = 0.00000    ; % P116: !=0: x elements for fast tensor product
i = i+1; param(i) = 0.00000    ; % P117: !=0: y elements for fast tensor product
i = i+1; param(i) = 0.00000    ; % P118: !=0: z elements for fast tensor product
