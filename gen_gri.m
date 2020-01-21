function gen_gri(gridname, nelx, nely,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Generate mesh for wing simulations (lower cut)                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - boundary/initial condition FLUENT data
dxbox = 2e-3; xbox = [-.2 .8]; nxbox = ceil(diff(xbox)/dxbox) + 1;
dybox = 2e-3; ybox = [-.2 .2]; nybox = ceil(diff(ybox)/dybox) + 1;

bcsource = 'FLUENT';
bcdata = 'naca0008_fluent-data.mat';


% - boundaries
% -- profile
xprcutup = .5;.3;
xprcutlw = .04;

% -- inflow
lin = .3;.25;.5;
xinup = -.05;-.00;-.05;-.00;
xinlw = -.05;-.00;-.05;-.00;

% -- outflows
doutup = .10;
doutlw = .08;

% -- fringe
wdfrup = .050; rsfrup = .015;
wdfrlw = .030; rsfrlw = .010;


% - grid

% -- exponent for weight function (how dx is related to the local curvature of the profile)
rexp = .55;.2;.25;1/2;


%% Load original profile coordinates
iaf.designation='0008';
iaf.n = n;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;
data = naca4gen(iaf); xpro =data.x'; ypro = flip(data.z)'; clear data

% tangential coordinate
dx = xpro(2:end) - xpro(1:end-1);
dy = ypro(2:end) - ypro(1:end-1);

ds = sqrt(dx.^2+dy.^2);
spro = cumsum([0 ds]); [~,ispro0] = min(abs(xpro)); spro = spro - spro(ispro0);
% profile spline
xprfun = csapi(spro,xpro);
yprfun = csapi(spro,ypro);


% curvature radius (spline differentiation)
d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
d1y = fnval(fnder(yprfun,1),spro); d2y = fnval(fnder(yprfun,2),spro);

rpro = abs(sqrt((d1x.^2 + d1y.^2).^3)./ ...
              (d1x.*d2y - d2x.*d1y)   );

rprfun = csapi(spro,rpro);

% tangential and normal directions
tpro = [  d1x ;
          d1y ]; tpro = tpro ./ ([1;1]*(sqrt(tpro(1,:).^2 + tpro(2,:).^2)));
npro = [ -d1y ;
          d1x ]; npro = npro ./ ([1;1]*(sqrt(npro(1,:).^2 + npro(2,:).^2)));
tprfun{1} = csapi(spro,tpro(1,:)); tprfun{2} = csapi(spro,tpro(2,:));
nprfun{1} = csapi(spro,npro(1,:)); nprfun{2} = csapi(spro,npro(2,:));

%% Cut and re-interpolate profile
% cut profile
icb = find(xpro < xprcutlw,1,'first');
ice = find(xpro < xprcutup,1,'last'); scut = spro([icb,ice]);

scut(1) = fzero(@(s) fnval(xprfun,s)-xprcutlw,scut(1));
scut(2) = fzero(@(s) fnval(xprfun,s)-xprcutup,scut(2));

% new spacing: weight function (wpr) for ds based on the profile curvature (1/rpr)
wprfun = fnint(csapi(spro,abs(1./rpro).^(rexp)));
wpro = fnval(wprfun,spro); wprinv = csapi(wpro,spro);

wpr = linspace(fnval(wprfun,scut(1)),fnval(wprfun,scut(2)),nelx+1);

% interpolation
spr = fnval(wprinv,wpr);
xpr = fnval(xprfun,spr);
ypr = fnval(yprfun,spr);
tpr =[fnval(tprfun{1},spr);
      fnval(tprfun{2},spr)];
npr =[fnval(nprfun{1},spr);
      fnval(nprfun{2},spr)];
rpr = fnval(rprfun,spr);

ipr = 1 + 0*xpr;

%% Mid points airfoil
spr_2 =  (spr(2:end) - spr(1:end-1))/2 + spr(1:end-1);
x_m = fnval(xprfun,spr_2)
y_m = fnval(yprfun,spr_2)
midpoints.x = x_m;
midpoints.y = y_m;
save('mid_points.mat', 'midpoints')
%% Fringe
Sfun = @(x) ( 1./(1 + exp(1./(x-1) + 1./x)) ).*(0 < x).*(x < 1) + (1 <= x);

% compute fringe (have a look at SIMSON manual for ref.)
frpr = Sfun((wdfrup- (spr(end)-spr))/rsfrup) + Sfun((wdfrlw - (spr-spr(1)))/rsfrlw);



%% Compute upper/lower boundary from FLUENT/Morino data (streamlines)
xic = linspace(xbox(1),xbox(2),nxbox);
yic = linspace(ybox(1),ybox(2),nybox);

[xxic,yyic] = meshgrid(xic,yic);

% interpolate data on a rectangular box (xf,yf)
bc = load(bcdata);

switch bcsource
    case 'Morino'
        [uuic,vvic,ppic] = morinoflowfield(bc.pan+bc.dpan,bc.cll+bc.dcll,...
                                                bc.alpha,bc.phi,xxic,yyic);

    case 'FLUENT'
        uuic = griddata(bc.xx,bc.yy,bc.uu,xxic,yyic,'cubic');
        vvic = griddata(bc.xx,bc.yy,bc.vv,xxic,yyic,'cubic');
        ppic = griddata(bc.xx,bc.yy,bc.pp,xxic,yyic,'cubic');
end

% streamline (upper boundaries)
x0up = xpr(end) + doutup*npr(1,end);
y0up = ypr(end) + doutup*npr(2,end);
line = stream2(xic,yic,-uuic,-vvic,x0up,y0up); xupo = line{1}(:,1)'; yupo = line{1}(:,2)';

% - tangential coordinate
dx = xupo(2:end) - xupo(1:end-1);
dy = yupo(2:end) - yupo(1:end-1);
ds = sqrt(dx.^2+dy.^2); supo = cumsum([0 ds]);

xupfun = csapi(supo,xupo); xupinv = csapi(xupo,supo);
yupfun = csapi(supo,yupo);


% streamline (lower boundaries)
x0lw = xpr(1) + doutlw*npr(1,1);
y0lw = ypr(1) + doutlw*npr(2,1);
line = stream2(xic,yic,-uuic,-vvic,x0lw,y0lw); xlwo = line{1}(:,1)'; ylwo = line{1}(:,2)';

% - tangential coordinate
dx = xlwo(2:end) - xlwo(1:end-1);
dy = ylwo(2:end) - ylwo(1:end-1);
ds = sqrt(dx.^2+dy.^2); slwo = cumsum([0 ds]);

xlwfun = csapi(slwo,xlwo); xlwinv = csapi(xlwo,slwo);
ylwfun = csapi(slwo,ylwo);



%% Connect and merge streamlines boundaries (freestream)

% cut upper boundary (by interpolation)
xup = linspace(xinup,xupo(1),ceil((xupo(1)-xinup)/dxbox)+1);
sup = fnval(xupinv,xup);
yup = fnval(yupfun,sup);
iup = 3 + 0*xup;


% cut lower boundary (by interpolation)
xlw = linspace(xinlw,xlwo(1),ceil((xlwo(1)-xinlw)/dxbox)+1);
slw = fnval(xlwinv,xlw);
ylw = fnval(ylwfun,slw);
ilw = 5 + 0*xlw;

% flip lower boundary for consistency
xlw = fliplr(xlw); ylw = fliplr(ylw); slw = fliplr(slw); slw = slw - slw(1);


% inflow (spline from upper to lower streamline)
xx = [         xlw(end-2:end) xup(1:3)           ];
yy = [         ylw(end-2:end) yup(1:3)           ];
ss = [slw(end-2:end)-slw(end) sup(1)-sup(1:3)+lin];

sin = linspace(0,lin,ceil(lin/dxbox)+1);
xin = interp1(ss,xx,sin,'spline');
yin = interp1(ss,yy,sin,'spline');
iin = 4 + 0*xin;


% merge streamlines and inflow (freestream)
xfso = [ xlw(1:end-1) xin(1:end) xup(2:end)];
yfso = [ ylw(1:end-1) yin(1:end) yup(2:end)];
ifso = [ ilw(1:end-1) iin(1:end) iup(2:end)];

% - tangential coordinate
dx = xfso(2:end) - xfso(1:end-1);
dy = yfso(2:end) - yfso(1:end-1);
ds = sqrt(dx.^2+dy.^2); sfso = cumsum([0 ds]);

% freestream spline
xfsfun = csapi(sfso,xfso);
yfsfun = csapi(sfso,yfso); 

% curvature radius (spline differentiation)
d1x = fnval(fnder(xfsfun,1),sfso); d2x = fnval(fnder(xfsfun,2),sfso);
d1y = fnval(fnder(yfsfun,1),sfso); d2y = fnval(fnder(yfsfun,2),sfso);

rfso = abs(sqrt((d1x.^2 + d1y.^2).^3)./ ...
              (d1x.*d2y - d2x.*d1y)   );

rfsfun = csapi(spro,rpro);

% tangential and normal directions
tfso = -[  d1x ;
           d1y ]; tfso = tfso ./ ([1;1]*(sqrt(tfso(1,:).^2 + tfso(2,:).^2)));
nfso = -[ -d1y ;
           d1x ]; nfso = nfso ./ ([1;1]*(sqrt(nfso(1,:).^2 + nfso(2,:).^2)));

tfsfun{1} = csapi(sfso,tfso(1,:)); tfsfun{2} = csapi(sfso,tfso(2,:));
nfsfun{1} = csapi(sfso,nfso(1,:)); nfsfun{2} = csapi(sfso,nfso(2,:));



%% Map the free-stream boundary via profile normals
% define tangential distance function
fun = @(s,t,x0) t' * ([fnval(xfsfun,s);
                       fnval(yfsfun,s)] - x0);

% map (by finding intersections)
sfs = zeros(size(spr));
for i = 2:length(spr)
    sfs(i) = fzero(@(s) fun(s,tpr(:,i),[xpr(i);
                                        ypr(i)]),sfs(i-1));
end
xfs = fnval(xfsfun,sfs);
yfs = fnval(yfsfun,sfs);

tfs =[fnval(tfsfun{1},sfs);
      fnval(tfsfun{2},sfs)];
nfs =[fnval(nfsfun{1},sfs);
      fnval(nfsfun{2},sfs)];
rfs = fnval(rfsfun,sfs);

ifs = interp1(sfso,ifso,sfs,'linear'); ifs = floor(ifs);



%% Build grid
nbd = 6; %number of boundaries

map = (1-cos(linspace(0,1,nely+1)*pi/2));

xxgr = zeros(nely+1,nelx+1);
yygr = zeros(nely+1,nelx+1);
iigr = zeros(nely+1,nelx+1,4); % 1=N, 2=E, 3=S, 4=W
frgr = zeros(nely+1,nelx+1);
%% curvature radius (spline differentiation)
%% curvature radius (spline differentiation)
d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);

for i = 1:nelx+1
    % wall-normal direction
    xxgr(:,i) = interp1([0 1],[xpr(i) xfs(i)],map,'linear');
    yygr(:,i) = interp1([0 1],[ypr(i) yfs(i)],map,'linear');

    % boundary id
    if i ~= nelx+1
        iigr(1  ,i,1) = ipr(i+1); % profile
        iigr(end,i,1) = ifs(i+1); % free-stream
    end
    if i ~= 1
        iigr(1  ,i,3) = ipr(i-1); % profile
        iigr(end,i,3) = ifs(i-1); % free-stream
    end
    if i == 1
        iigr(2:end  ,i,2) = 6; % lower outflow boundary id
        iigr(1:end-1,i,4) = 6; % lower outflow boundary id
    elseif i == nelx+1
%% curvature radius (spline differentiation)
d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
        iigr(2:end  ,i,2) = 2; % upper outflow boundary id
        iigr(1:end-1,i,4) = 2; % upper outflow boundary id
    end

    % curvature radius (profile only)
    frgr(:,i) = frpr(i);
end



%% Initial values and boundary conditions

switch bcsource
    case 'Morino'
        [uugr,vvgr,ppgr] = morinoflowfield(bc.pan+bc.dpan,bc.cll+bc.dcll,bc.alpha,bc.phi,xxgr,yygr);

        % fix boundary layer
        bl = load(bldata);
        dstfun = @(m) interp1(bl.m,bl.dst,m,'spline');
        ublfun = @(eta,m) interp2(bl.m,bl.yy,squeeze(bl.ff(2,:,:)),m,eta,'spline',1);
        etamax = min([max(bl.yy) 20])

        % wall normal direction
        n = [xxgr(end,:) - xxgr(1,:);
             yygr(end,:) - yygr(1,:)]; n = n./([1 1]' * sqrt(n(1,:).^2 + n(2,:).^2));

        % boundary layer characteristics
        mbl = interp1(bc.scll,bc.m,spr,'spline');
        deltabl = interp1(bc.scll,bc.delta,spr,'spline');

        for i = 1:nelx+1

            % boundary layer coordinate
            etabl = n(:,i)' * [ xxgr(:,i)' - xxgr(1,i) ;
                                yygr(:,i)' - yygr(1,i) ] ./ deltabl(i);
            jj = etabl <= etamax;

            % regularize velocity (close to the wall)
            method = 'linear';
            uugr(jj,i) = interp1(etabl(~jj),uugr(~jj,i),etabl(jj),method,'extrap');
            vvgr(jj,i) = interp1(etabl(~jj),vvgr(~jj,i),etabl(jj),method,'extrap');
            ppgr(jj,i) = 1 - (uugr(jj,i).^2 + vvgr(jj,i).^2)/2;

            % map boundary layer to grid
            uugr(jj,i) = ublfun(etabl(jj),mbl(i)) .* uugr(jj,i);
            vvgr(jj,i) = ublfun(etabl(jj),mbl(i)) .* vvgr(jj,i);

        end

    case 'FLUENT'
        uugr = griddata(bc.xx,bc.yy,bc.uu,xxgr,yygr,'cubic');
        vvgr = griddata(bc.xx,bc.yy,bc.vv,xxgr,yygr,'cubic');
        ppgr = griddata(bc.xx,bc.yy,bc.pp,xxgr,yygr,'cubic');

        %Pgr = ConstructPolyInterpolant2D(bc.xx,bc.yy,xxgr,yygr,3,-64);

        %uugr = reshape(Pgr*bc.uu(:),size(xxgr));
        %vvgr = reshape(Pgr*bc.vv(:),size(xxgr));
        %ppgr = reshape(Pgr*bc.pp(:),size(xxgr));
end

% fix no-slip on the profile
%bd = sum(iigr==1,3)~=0; 
%uugr(bd) = 0;
%vvgr(bd) = 0;



%% Check mass-flux
% - profile
prflux = trapz(spr,uugr(1,:).*npr(1,:) + vvgr(1,:).*npr(2,:))
% - upper outflow
luo = sqrt((xpr(end)-xfs(end))^2 + (ypr(end)-yfs(end))^2);
uoflux = trapz(map*luo,-uugr(:,end).*tpr(1,end) - vvgr(:,end).*tpr(2,end))
% - free-stream
ii = ifs==3; upflux = trapz(sfs(ii),uugr(end,ii).*nfs(1,ii) + vvgr(end,ii).*nfs(2,ii))
ii = ifs==4; influx = trapz(sfs(ii),uugr(end,ii).*nfs(1,ii) + vvgr(end,ii).*nfs(2,ii))
ii = ifs==5; lwflux = trapz(sfs(ii),uugr(end,ii).*nfs(1,ii) + vvgr(end,ii).*nfs(2,ii))
fsflux = trapz(sfs,uugr(end,:).*nfs(1,:) + vvgr(end,:).*nfs(2,:))
% - lower outflow
llo = sqrt((xpr(1)-xfs(1))^2 + (ypr(1)-yfs(1))^2);
loflux = trapz(map*llo, uugr(:,1).*tpr(1,1) + vvgr(:,1).*tpr(2,1))
% - TOTAL
flux = prflux + fsflux + loflux + uoflux



%% Grid-quality check
xxel = zeros(nely,nelx); yyel = zeros(nely,nelx);
dtel = zeros(nely,nelx); dnel = zeros(nely,nelx);
arel = zeros(nely,nelx); skel = zeros(nely,nelx);

for i = 1:nelx
    for j = 1:nely
        v1 = [xxgr(j  ,i  ); yygr(j  ,i  )];
        v2 = [xxgr(j  ,i+1); yygr(j  ,i+1)];
        v3 = [xxgr(j+1,i+1); yygr(j+1,i+1)];
        v4 = [xxgr(j+1,i  ); yygr(j+1,i  )];

        l1 = v2-v1; ll1 = norm(l1);
        l2 = v3-v2; ll2 = norm(l2);
        l3 = v4-v3; ll3 = norm(l3);
        l4 = v1-v4; ll4 = norm(l4);

        th1 = acos((l1(1)*l4(1) + l1(2)*l4(2))/(ll1*ll4));
        th2 = acos((l2(1)*l1(1) + l2(2)*l1(2))/(ll2*ll1));
        th3 = acos((l3(1)*l2(1) + l3(2)*l2(2))/(ll3*ll2));
        th4 = acos((l4(1)*l3(1) + l4(2)*l3(2))/(ll4*ll3));

        d1 = v3 - v1;
        d2 = v4 - v2;

% element grid (element center of mass)
        xxel(j,i) = (v1(1) + v2(1) + v3(1) + v4(1))/4;
        yyel(j,i) = (v1(2) + v2(2) + v3(2) + v4(2))/4;

% aspect ratio
        arel(j,i) = max([ll1 ll2 ll3 ll4])/min([ll1 ll2 ll3 ll4]);

% skewness
        skel(j,i) = max(abs([th1 th2 th3 th4] - pi/2))/(pi/2);

% resolution
        dnel(j,i) = max([ll2 ll4]);
        dtel(j,i) = max([ll1 ll3]);


    end
end



%% Save grid
write_grid(gridname,xxgr,yygr,iigr,uugr,vvgr,ppgr,frgr);



%% %% Plot grid
%% figure(1); clf; hold on; set (1,'Units','normalized','Position',[.5 0 .5 1]);

%% % ppic(abs(ppic)>3) = NaN;
%% % hh = surf(xxic,yyic,uuic,'EdgeColor','none'); hc = colorbar('NO'); xlabel(hc,'u')
%% hh = surf(xxic,yyic,ppic,'EdgeColor','none'); hc = colorbar('NO'); xlabel(hc,'p')
%% os = max(max(get(hh,'ZData')));

%% plot3(xpro,ypro,os+0*xpro,'--w','LineWidth',2.0); % profile
%% plot3(xlwo,ylwo,os+0*xlwo,'--w','LineWidth',2.0); % lower streamline
%% plot3(xupo,yupo,os+0*xupo,'--w','LineWidth',2.0); % upper streamline

%% hgr = mesh(xxgr,yygr,os+0*xxgr,'EdgeColor',[1 1 1],'FaceColor','none'); % grid

%% bd = sum(iigr==1,3)~=0; hbd(1) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-r','LineWidth',1.5);
%% bd = sum(iigr==2,3)~=0; hbd(2) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-g','LineWidth',1.5);
%% bd = sum(iigr==3,3)~=0; hbd(3) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-b','LineWidth',1.5);
%% bd = sum(iigr==4,3)~=0; hbd(4) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-c','LineWidth',1.5);
%% bd = sum(iigr==5,3)~=0; hbd(5) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-m','LineWidth',1.5);
%% bd = sum(iigr==6,3)~=0; hbd(6) = plot3(xxgr(bd),yygr(bd),os+0*xxgr(bd),'-y','LineWidth',1.5);

%% legend(hbd,'1: profile',...
%%            '2: upper outflow',...
%%            '3: upper streamline',...
%%            '4: inflow',...
%%            '5: lower stremline',...
%%            '6: lower outflow','Location','SO','Orientation','Horizontal');

%% hold off; view(2); axis equal; axis([xbox ybox]);
%% xlabel('x/c'); ylabel('y/c'); title('Grid'); drawnow



%% figure(2); clf; set (2,'Units','normalized','Position',[0 0 .5 1]);

%% subplot(4,1,1); hh = surf(xxgr,yygr,uugr,'EdgeColor','none'); hc = colorbar('EO'); xlabel(hc,'u'); hold on
%%                 if bcsource == 'Morino'
%%                 plot3(xxgr(1,:)+etamax*deltabl.*n(1,:),yygr(1,:)+etamax*deltabl.*n(2,:),xxgr(1,:)*0+max(max(uugr)),'w');
%%                 end;
%%                 hold off
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title('RANS x-velocity')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxgr(inmax,itmax),yygr(inmax,itmax),os,'^r',...
%%                            xxgr(inmin,itmin),yygr(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,2); hh = surf(xxgr,yygr,vvgr,'EdgeColor','none'); hc = colorbar('EO'); xlabel(hc,'v'); hold on
%%                 if bcsource == 'Morino'
%%                 plot3(xxgr(1,:)+etamax*deltabl.*n(1,:),yygr(1,:)+etamax*deltabl.*n(2,:),xxgr(1,:)*0+max(max(vvgr)),'w');
%%                 end;
%%                 hold off
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title('RANS y-velocity')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxgr(inmax,itmax),yygr(inmax,itmax),os,'^r',...
%%                            xxgr(inmin,itmin),yygr(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,3); hh = surf(xxgr,yygr,ppgr,'EdgeColor','none'); hc = colorbar('EO'); xlabel(hc,'p'); hold on
%%                 if bcsource == 'Morino'
%%                 plot3(xxgr(1,:)+etamax*deltabl.*n(1,:),yygr(1,:)+etamax*deltabl.*n(2,:),xxgr(1,:)*0+max(max(ppgr)),'w');
%%                 end;
%%                 hold off
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title('RANS pressure')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxgr(inmax,itmax),yygr(inmax,itmax),os,'^r',...
%%                            xxgr(inmin,itmin),yygr(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,4); hh = surf(xxgr,yygr,frgr,'EdgeColor','none'); hc = colorbar('EO'); xlabel(hc,'f'); hold on
%%                 if bcsource == 'Morino'
%%                 plot3(xxgr(1,:)+etamax*deltabl.*n(1,:),yygr(1,:)+etamax*deltabl.*n(2,:),xxgr(1,:)*0+max(max(frgr)),'w');
%%                 end;
%%                 hold off
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title('Fringe-forcing mask')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxgr(inmax,itmax),yygr(inmax,itmax),os,'^r',...
%%                            xxgr(inmin,itmin),yygr(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off



%% figure(3); clf; set (3,'Units','normalized','Position',[.5 0 .5 1]);

%% subplot(4,1,1); hh = surf(xxel,yyel,arel,'EdgeColor','none');
%%                 hc = colorbar('EO'); xlabel(hc,'|l_i|_{max} / |l_i|_{min}')
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title ('Aspect Ratio')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxel(inmax,itmax),yyel(inmax,itmax),os,'^r',...
%%                            xxel(inmin,itmin),yyel(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,2); hh = surf(xxel,yyel,skel,'EdgeColor','none');
%%                 hc = colorbar('EO'); xlabel(hc,'|\theta_i - \pi/2|_{max} / \pi/2')
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title ('Skewness')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxel(inmax,itmax),yyel(inmax,itmax),os,'^r',...
%%                            xxel(inmin,itmin),yyel(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,3); hh = surf(xxel,yyel,dtel,'EdgeColor','none');
%%                 hc = colorbar('EO'); xlabel(hc,'\Deltat')
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title ('Resolution (wall-wise)')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxel(inmax,itmax),yyel(inmax,itmax),os,'^r',...
%%                            xxel(inmin,itmin),yyel(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off

%% subplot(4,1,4); hh = surf(xxel,yyel,dnel,'EdgeColor','none');
%%                 hc = colorbar('EO'); xlabel(hc,'\Deltan')
%%                 view(2); axis image; grid on
%%                 xlabel('x/c'); ylabel('y/c'); title ('Resolution (wall-normal)')

%%                 ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
%%                 [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
%%                 [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
%%                 hm = plot3(xxel(inmax,itmax),yyel(inmax,itmax),os,'^r',...
%%                            xxel(inmin,itmin),yyel(inmin,itmin),os,'vb','MarkerFaceColor','w');
%%                 legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%%                           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off


% figure(4); clf; set (4,'Units','normalized','Position',[0. 0. 1. 1.]);
% hh = surf(xxgr,yygr,uugr,'EdgeColor','none'); hc = colorbar('EO'); xlabel(hc,'u'); hold on
% plot3(xxgr(1,:)+etamax*deltabl.*n(1,:),yygr(1,:)+etamax*deltabl.*n(2,:),xxgr(1,:)*0+max(max(uugr)),'w'); hold off
% view(2); axis image; grid on
% xlabel('x/c'); ylabel('y/c'); title('RANS x-velocity')
%
% ax = axis; hold on; data = get(hh,'ZData'); os = max(max(data));
% [dum,inmax] = max(data); [dum,itmax] = max(dum); inmax = inmax(itmax);
% [dum,inmin] = min(data); [dum,itmin] = min(dum); inmin = inmin(itmin);
% hm = plot3(xxgr(inmax,itmax),yygr(inmax,itmax),os,'^r',...
%            xxgr(inmin,itmin),yygr(inmin,itmin),os,'vb','MarkerFaceColor','w');
% legend(hm,sprintf('max: %.2e',data(inmax,itmax)),...
%           sprintf('min:  %.2e',data(inmin,itmin)),'Location','NE'); hold off
