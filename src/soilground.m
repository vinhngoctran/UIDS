% OFM-Urban, original author Vinh Ngoc Tran (2021).
% Soil/vegetation processes follow the BROOK90 / Hammel & Kennel SVAT model.
% Hydraulic conventions and repair notes: docs/HYDROLOGY.md.
function [psim,gwat,vrfli,infli,byfli,dsfli,ntfli,trani,swat,gwfl,seep,diagnostics] = ...
    soilground(p,psim,ptr,gwat,slvp,slfl,deltaT)
%SOILGROUND Layered BROOK90-style soil/vegetation water balance.
% Input demands ptr/slvp and infiltration slfl are m/s; timestep is seconds.
% State: matric potentials psim [kPa, <=0] and groundwater gwat [mm].
% Output layer rates and gwfl/seep are timestep means [mm/day].
% imodel=0: modified Clapp-Hornberger (legacy default); imodel=1: Mualem-VG.
% All 18 original process helpers remain, with consistent hydraulic curves.
validateattributes(deltaT,{'numeric'},{'scalar','positive','finite'});
for value={ptr,slvp,slfl,gwat}
    validateattributes(value{1},{'numeric'},{'scalar','real','finite','nonnegative'});
end
n = p.nlayer;
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
imodel = field(p,'imodel',0);
validateattributes(field(p,'max_step',300),{'numeric'},{'scalar','real','finite','positive'});
validateattributes(field(p,'dpsimx',0.001),{'numeric'},{'scalar','real','finite','nonnegative'});
thick = layer(p.thick,n); stone = layer(p.stonef,n);
psim = layer(psim,n);
validateattributes(psim,{'numeric'},{'real','finite'});
if any(~isfinite(psim) | psim > 0), error('OFM:SoilPotential','Matric potential must be finite and nonpositive [kPa].'); end
[psig,capacity,~,~,~,~,wet,water,ksat,par] = soilpar(n,p.par,thick,stone,psim,p.psicr,imodel);
minimum = zeros(1,n);
if imodel == 1, minimum = thick.*(1-stone).*par(10,:); end
initial = sum(water)+gwat;
startWater = water;
days = deltaT/86400;
inputRate = slfl*1000*86400;
potentialT = ptr*1000*86400; potentialE = slvp*1000*86400;
infrac = infpar(p.infexp,p.ilayer,n,thick);
relden = rootgrowth(layer(p.frelden,n),layer(p.tini,n),p.age,p.rgroper,p.inirdep,p.inirlen,n);
[~,~,~,rtlen,rplant] = canopy(p.snow,p.snoden,p.mxrtln,p.mxkpl,p.densef,p.height,p.lai,p.sai);
[rxylem,rrooti,alpha] = plntres(n,thick,stone,rtlen,relden,p.rtrad,rplant,p.fxylem);
vrfli=zeros(1,n); infli=vrfli; byfli=vrfli; dsfli=vrfli; trani=vrfli;
gwDepth=0; seepDepth=0; evapDepth=0;
remaining=days; substeps=0;
while remaining > 0
    substeps=substeps+1;
    if substeps > 100000, error('OFM:SoilSubsteps','Too many soil substeps. Reduce DELTA_T or check soil parameters.'); end
    [psiti,~,~,kk] = soilvar(n,par,psig,psim,wet,water,imodel);
    if ~isreal(kk) || any(~isfinite(kk) | kk<0)
        error('OFM:SoilConductivity','Hydraulic parameters produced invalid conductivity.');
    end
    q=zeros(1,n);
    for i=1:n-1
        if abs(psiti(i)-psiti(i+1)) > field(p,'dpsimx',0.001)
            q(i)=vert(kk(i),kk(i+1),ksat(i),ksat(i+1),thick(i),thick(i+1), ...
                psiti(i),psiti(i+1),stone(i),stone(i+1),0.00981);
        end
    end
    q(n)=max(0,p.drain)*kk(n)*(1-stone(n));
    down=dslop(p.dslope,p.length,thick,stone,psim,0.00981,kk);
    [~,uptake]=tbylayer(potentialT,p.dispc,alpha,kk,rrooti,rxylem,psiti,n,-abs(p.psicr));
    outgoing=down+uptake;
    outgoing(1)=outgoing(1)+potentialE;
    for i=1:n-1
        if q(i)>=0, outgoing(i)=outgoing(i)+q(i); else, outgoing(i+1)=outgoing(i+1)-q(i); end
    end
    outgoing(n)=outgoing(n)+q(n);
    dt=min(remaining,field(p,'max_step',300)/86400);
    mobile=water-minimum;
    draining=mobile>1e-9 & outgoing>0;
    if any(draining), dt=min(dt,0.25*min(mobile(draining)./outgoing(draining))); end
    if ~isfinite(dt) || dt<=0 || remaining-dt==remaining, error('OFM:SoilSubsteps','Soil timestep underflow.'); end
    byfrac=byflfr(n,wet,par,p.qffc,p.qfpar);
    if field(p,'bypar',1)==0, byfrac(:)=0; end
    % Apply infiltration once. Excess above capacity remains bypass runoff.
    gross=inputRate*infrac*dt;
    bypass=gross.*byfrac;
    admitted=min(gross-bypass,max(0,capacity-water));
    bypass=gross-admitted;
    water=water+admitted;
    infli=infli+admitted; byfli=byfli+bypass;
    % Evaporation and transpiration remove only available mobile soil water.
    actualE=min(potentialE*dt,max(0,water(1)-minimum(1)));
    water(1)=water(1)-actualE; evapDepth=evapDepth+actualE;
    actualT=min(uptake*dt,max(0,water-minimum));
    water=water-actualT; trani=trani+actualT;
    actualD=min(down*dt,max(0,water-minimum));
    water=water-actualD; dsfli=dsfli+actualD;
    % Paired layer transfers preserve volume and both storage bounds.
    for i=1:n-1
        if q(i)>=0, donor=i; receiver=i+1; else, donor=i+1; receiver=i; end
        v=min([abs(q(i))*dt,max(0,water(donor)-minimum(donor)),max(0,capacity(receiver)-water(receiver))]);
        water(donor)=water(donor)-v; water(receiver)=water(receiver)+v;
        vrfli(i)=vrfli(i)+sign(q(i))*v;
    end
    drainageDepth=min(q(n)*dt,max(0,water(n)-minimum(n)));
    water(n)=water(n)-drainageDepth;
    vrfli(n)=vrfli(n)+drainageDepth;
    [gflow,deep]=gwater(gwat,p.gsc,p.gsp,dt,drainageDepth/dt);
    gwat=max(0,gwat+drainageDepth-(gflow+deep)*dt);
    gwDepth=gwDepth+gflow*dt; seepDepth=seepDepth+deep*dt;
    wet=storage_wetness(water,capacity,minimum,imodel);
    for i=1:n, psim(i)=fpsim(wet(i),par(:,i),imodel); end
    remaining=max(0,remaining-dt);
end
swat=sum(water);
ntfli=(water-startWater)/days;
vrfli=vrfli/days; infli=infli/days; byfli=byfli/days;
dsfli=dsfli/days; trani=trani/days;
gwfl=gwDepth/days; seep=seepDepth/days;
outputs=evapDepth+sum(trani)*days+sum(dsfli+byfli)*days+gwDepth+seepDepth;
balance=initial+inputRate*days-swat-gwat-outputs;
if ~isreal([psim gwat balance]) || any(~isfinite([psim gwat balance outputs]))
    error('OFM:SoilState','Soil integration produced a nonfinite or complex state.');
end
if abs(balance)>1e-8*max(1,initial+inputRate*days)
    error('OFM:SoilBalance','Soil water balance residual is %.12g mm.',balance);
end
diagnostics=struct('SoilEvaporation',evapDepth/days,'Transpiration',sum(trani), ...
    'BypassRunoff',sum(byfli),'DownslopeRunoff',sum(dsfli),'BalanceResidual',balance, ...
    'InitialStorage',initial,'LayerStorage',water,'MinimumStorage',minimum, ...
    'Capacity',capacity,'Substeps',substeps);
end

function value=field(p,name,default)
if isfield(p,name), value=p.(name); else, value=default; end
end

function values=layer(values,n)
values=double(values(:).');
if isscalar(values), values=repmat(values,1,n); end
if numel(values)~=n, error('OFM:SoilLayers','A layer parameter/state has the wrong length.'); end
end

function wet=storage_wetness(water,capacity,minimum,imodel)
if imodel==0, wet=water./capacity; else, wet=(water-minimum)./(capacity-minimum); end
wet=max(0,min(1,wet));
end

function [psig,swatmx,wetf,wetc,chm,chn,wetnes,swati,ksat,par] = soilpar(nlayer,par,thick,stonef,psim,psicr,imodel)
if nargin<7, imodel=0; end
if ~ismember(imodel,[0 1]) || ~isequal(size(par),[10 nlayer])
    error('OFM:SoilParameters','Expected a 10-by-NLAYER parameter matrix and imodel 0 or 1.');
end
if ~isreal(par) || any(~isfinite([par(:);thick(:);stonef(:)])) || any(thick<=0) || any(stonef<0 | stonef>=1)
    error('OFM:SoilParameters','Invalid layer thickness, stone fraction or hydraulic parameters.');
end
if any(par(1,:)<=0 | par(1,:)>1), error('OFM:SoilParameters','Saturated water contents must lie in (0,1].'); end
psig=-0.00981*(cumsum(thick)-thick/2);
swatmx=thick.*par(1,:).*(1-stonef);
wetnes=zeros(1,nlayer); swati=wetnes; wetf=wetnes; wetc=wetnes; ksat=wetnes; chm=wetnes; chn=wetnes;
for i=1:nlayer
    if imodel==0
        if par(2,i)<=0 || par(2,i)>=par(1,i) || par(3,i)<=0 || par(4,i)==0 || par(9,i)<=0
            error('OFM:SoilParameters','Invalid Clapp-Hornberger field capacity, conductivity, potential or exponent.');
        end
        par(4,i)=-abs(par(4,i));
        wetf(i)=par(2,i)/par(1,i); wi=par(10,i); b=par(9,i);
        if wi<=wetf(i) || wi>=1, error('OFM:SoilParameters','wetinf must exceed field-capacity saturation and be below 1.'); end
        par(5,i)=wetf(i);
        psiinf=par(4,i)*(wi/wetf(i))^(-b);
        chm(i)=-psiinf/(1-wi)^2+b*psiinf/(wi*(1-wi));
        if chm(i)<=0, error('OFM:SoilParameters','Near-saturation parameters do not form a monotonic retention curve.'); end
        chn(i)=2*wi-1+psiinf*b/(chm(i)*wi);
        par(7,i)=chm(i); par(8,i)=chn(i);
        ksat(i)=par(3,i)*(1/wetf(i))^(2*b+3); par(6,i)=ksat(i);
    else
        if par(10,i)<0 || par(10,i)>=par(1,i) || par(6,i)<=0 || par(7,i)<=0 || par(8,i)<=1
            error('OFM:SoilParameters','Invalid Mualem-van Genuchten residual content, Ksat, alpha or n.');
        end
        ksat(i)=par(6,i);
        if par(3,i)<=0 || par(3,i)>=ksat(i)
            error('OFM:SoilParameters','VG field-capacity conductivity must lie between 0 and Ksat.');
        end
        lo=0; hi=1;
        for k=1:60
            mid=(lo+hi)/2;
            if fk(mid,par(:,i),imodel)<par(3,i), lo=mid; else, hi=mid; end
        end
        wetf(i)=(lo+hi)/2; par(5,i)=wetf(i);
        par(4,i)=fpsim(wetf(i),par(:,i),imodel);
        par(2,i)=par(10,i)+wetf(i)*(par(1,i)-par(10,i));
    end
    wetnes(i)=wet_from_potential(psim(i),par(:,i),imodel);
    wetc(i)=wet_from_potential(-abs(psicr)*1000,par(:,i),imodel);
    if imodel==0
        swati(i)=wetnes(i)*swatmx(i);
    else
        theta=par(10,i)+wetnes(i)*(par(1,i)-par(10,i));
        swati(i)=thick(i)*(1-stonef(i))*theta;
    end
end
if any(~isfinite(ksat)) || any(ksat<=0)
    error('OFM:SoilParameters','Derived saturated conductivities must be finite and positive.');
end
end

function wet=wet_from_potential(psi,p,imodel)
if psi>=0, wet=1; return; end
if psi<=-realmax, wet=0; return; end
if imodel==0
    wet=exp(log(p(5))-(log(-psi)-log(-p(4)))/p(9));
    if wet>p(10)
        lo=p(10); hi=1;
        for k=1:60
            mid=(lo+hi)/2;
            if fpsim(mid,p,imodel)<psi, lo=mid; else, hi=mid; end
        end
        wet=(lo+hi)/2;
    end
else
    m=1-1/p(8);
    exponent=p(8)*(log(p(7))+log(-psi)-log(9.81));
    softplus=max(0,exponent)+log1p(exp(-abs(exponent)));
    wet=exp(-m*softplus);
end
wet=max(0,min(1,wet));
end

function [psiti,theta,swat,kk]=soilvar(n,par,psig,psim,wetnes,swati,imodel)
if nargin<7, imodel=0; end
psiti=psim+psig; theta=zeros(1,n); kk=theta;
for i=1:n
    if imodel==0, theta(i)=wetnes(i)*par(1,i);
    else, theta(i)=par(10,i)+wetnes(i)*(par(1,i)-par(10,i)); end
    kk(i)=fk(wetnes(i),par(:,i),imodel);
end
swat=sum(swati);
end

function [vv,byfli,infli,ntfli]=inflow(n,dt,infrac,byfrac,slfl,vrfli,dsfli,trani,slvp,capacity,water)
% Reference flux assembly retained for older integrations; all fluxes mm/day.
infli=slfl*infrac.*(1-byfrac); byfli=slfl*infrac-infli; vv=vrfli;
for i=n:-1:1
    incoming=0; if i>1, incoming=vv(i-1); end
    demand=dsfli(i)+trani(i); if i==1, demand=demand+slvp; end
    maxInput=max(0,(capacity(i)-water(i))/dt+vv(i)+demand-incoming);
    admitted=min(infli(i),maxInput); byfli(i)=byfli(i)+infli(i)-admitted; infli(i)=admitted;
end
ntfli=infli+[0 vv(1:end-1)]-vv-dsfli-trani; ntfli(1)=ntfli(1)-slvp;
end

function byfrac=byflfr(n,wetnes,par,qffc,qfpar)
qffc=layer(qffc,n); qfpar=layer(qfpar,n);
if any(qffc<0 | qffc>1 | qfpar<=0), error('OFM:SoilBypass','Invalid bypass parameters.'); end
byfrac=zeros(1,n);
for i=1:n
    if qffc(i)==0, continue; end
    exponent=1-(wetnes(i)-par(5,i))/(qfpar(i)*max(1e-12,1-par(5,i)));
    byfrac(i)=exp(min(0,log(qffc(i))*exponent));
    if wetnes(i)>0.99, byfrac(i)=1; end
end
end

function [height,lai,sai,rtlen,rplant]=canopy(snow,snoden,mxrtln,mxkpl,densef,height,lai,sai)
if snow>0 && snoden<=0, error('OFM:SnowDensity','Positive snow requires positive snow density.'); end
snowDepth=0; if snow>0, snowDepth=0.001*snow/snoden; end
original=max(0.01,height); ratio=max(0,original-snowDepth)/original;
height=max(0.01,original-snowDepth); lai=max(1e-5,ratio*densef*lai); sai=max(0,densef*sai);
rtlen=max(0,densef*mxrtln); rplant=1/max(1e-20,densef*mxkpl);
end

function dsfli=dslop(slope,length,thick,stone,psim,rhowg,kk)
if length<=0, error('OFM:SlopeLength','Slope length must be positive.'); end
ll=1000*length;
grad=rhowg*sin(slope)+2*(psim/ll)*cos(slope);
ratio=thick.*(1-stone)*max(0,cos(slope))/ll;
dsfli=kk.*ratio.*max(0,grad/rhowg);
end

function value=fk(wet,p,imodel)
wet=max(0,min(1,wet));
if wet==0, value=0; return; end
if imodel==0
    value=p(3)*(wet/p(5))^(2*p(9)+3);
else
    m=1-1/p(8);
    logWet=log(wet); powered=exp(logWet/m);
    if powered==0
        logTerm=log(m)+logWet/m;
    else
        logTerm=log(-expm1(m*log1p(-powered)));
    end
    value=p(6)*exp(p(9)*logWet+2*logTerm);
end
end

function psi=fpsim(wet,p,imodel)
if nargin<3, imodel=0; end
wet=max(0,min(1,wet));
if wet==1, psi=0; return; end
if wet==0, psi=-realmax; return; end
if imodel==0
    if wet<p(10), psi=-exp(min(log(realmax),log(-p(4))-p(9)*(log(wet)-log(p(5)))));
    else, psi=p(7)*(wet-p(8))*(wet-1); end
else
    m=1-1/p(8);
    exponent=-log(wet)/m;
    logHead=(exponent+log1p(-exp(-exponent)))/p(8)-log(p(7))+log(9.81);
    psi=-exp(min(log(realmax),logHead));
end
psi=max(-realmax,psi);
end

function [gwfl,seep]=gwater(gwat,gsc,gsp,dt,recharge)
if gsc<0 || gsp<0 || gsp>1, error('OFM:Groundwater','Invalid groundwater recession or seepage fraction.'); end
if gsc<1e-8
    total=recharge;
else
    fraction=-expm1(-gsc*dt);
    final=gwat*(1-fraction)+recharge/gsc*fraction;
    total=max(0,(gwat+recharge*dt-final)/dt);
end
seep=gsp*total; gwfl=(1-gsp)*total;
end

function fraction=infpar(exponent,ilayer,n,thick)
if ilayer<1 || ilayer>n || ilayer~=round(ilayer), error('OFM:SoilLayers','ilayer must lie in 1:NLAYER.'); end
fraction=zeros(1,n);
if exponent<=0, fraction(1)=1; return; end
edges=[0 cumsum(thick(1:ilayer))]/sum(thick(1:ilayer));
fraction(1:ilayer)=diff(edges.^exponent);
end

function value=interp(npairs,funct,xvalue,xx,yy) %#ok<INUSD>
x=funct(1:2:2*npairs); y=funct(2:2:2*npairs);
if any(diff(x)<=0), error('OFM:Interpolation','Interpolation coordinates must increase.'); end
if npairs==1, value=y(1); else, value=interp1(x,y,max(x(1),min(x(end),xvalue)),'linear'); end
end

function dtnew=iter(n,dt,dpsidw,ntfli,capacity,psiti,dswmax,dpsimx,thick,wet,par,imodel,trani,slvp,mpar,ml,dtmin,pr) %#ok<INUSD>
% Retained reference timestep limiter. No calls with incompatible signatures.
dtnew=dt;
active=abs(ntfli)>0;
if any(active), dtnew=min(dtnew,min(0.01*dswmax*capacity(active)./abs(ntfli(active)))); end
if imodel==0
    potentialRate=ntfli.*dpsidw./capacity;
else
    potentialRate=ntfli.*dpsidw./(capacity.*(1-par(10,:)./par(1,:)));
end
for i=1:n
    if ntfli(i)<0
        if imodel==0, mobile=wet(i)*capacity(i);
        else, mobile=wet(i)*capacity(i)*(1-par(10,i)/par(1,i)); end
        if mobile>0, dtnew=min(dtnew,0.5*mobile/(-ntfli(i))); end
    end
    if i<n
        gradient=psiti(i)-psiti(i+1);
        change=potentialRate(i)-potentialRate(i+1);
        predicted=gradient+dt*change;
        if abs(predicted)>dpsimx && gradient*predicted<0
            dtnew=min(dtnew,-gradient/change);
        end
    end
end
dtnew=max(dtmin,dtnew);
end

function [rxylem,rrooti,alpha]=plntres(n,thick,stone,rtlen,relden,rtrad,rplant,fxylem)
if fxylem<0 || fxylem>1 || rtrad<=0, error('OFM:Roots','Invalid root radius or xylem fraction.'); end
rxylem=fxylem*rplant; effective=thick.*(1-stone); total=sum(relden.*effective);
rrooti=ones(1,n)*1e20; alpha=rrooti;
if total<=0 || rtlen<0.1, return; end
for i=1:n
    if relden(i)<1e-5, continue; end
    fraction=relden(i)*effective(i)/total;
    rrooti(i)=(rplant-rxylem)/fraction;
    density=fraction*0.001*rtlen/effective(i);
    delta=pi*rtrad^2*density;
    if abs(delta-1)<1e-6, geometric=0; else, geometric=delta-3-2*log(delta)/(1-delta); end
    alpha(i)=max(0,geometric/(8*pi*density)*0.001*0.00981/effective(i));
end
end

function density=rootgrowth(final,start,age,period,depth,length,n)
final=layer(final,n); start=layer(start,n); density=final;
if period<=0, return; end
if depth<=0 || length<=0, error('OFM:Roots','Initial root depth and length must be positive during growth.'); end
initial=length/depth;
for i=1:n
    if age<start(i), density(i)=0;
    elseif age<start(i)+period
        fraction=(age-start(i))/period;
        density(i)=initial*(final(i)/initial)^fraction;
    end
end
end

function [maximum,fieldStorage]=srfpar(qlayer,par,thick,stone,capacity)
indices=1:qlayer;
maximum=sum(capacity(indices)); fieldStorage=sum(par(2,indices).*thick(indices).*(1-stone(indices)));
end

function rate=swge(aa,asubs,vpd,raa,ras,rss,delta,arate)
% Shuttleworth-Wallace ground evaporation, flux input/output in mm/day.
gamma=0.066; cprho=1220; conversion=86400/2.45e6;
lec=arate/conversion;
rs=(delta+gamma)*ras+gamma*rss; ra=(delta+gamma)*raa;
if rs+ra<=0, error('OFM:Resistance','Evaporation resistance must be positive.'); end
le=rs/(rs+ra)*lec+(cprho*vpd+delta*ras*asubs+delta*raa*aa)/(rs+ra);
rate=max(0,conversion*(le-lec));
end

function [atr,atrani]=tbylayer(ptr,dispc,alpha,kk,rrooti,rxylem,psiti,n,psicr)
atr=0; atrani=zeros(1,n);
if ptr<=0, return; end
flag=rrooti>1e15 | psiti/1000<=psicr | kk<=0;
for repeat=1:n+1
    valid=~flag;
    if ~any(valid), return; end
    resistance=max(1e-12,rrooti+alpha./max(kk,realmin));
    total=1/sum(1./resistance(valid));
    potential=total*sum(psiti(valid)./resistance(valid));
    supply=(potential/1000-psicr-0.00981*dispc)/(total+rxylem);
    r=(2/pi)*max(0,supply)/ptr;
    if r<1, atr=ptr*(1+r*acos(r)-sqrt(max(0,1-r^2))); else, atr=ptr; end
    atrani(:)=0;
    atrani(valid)=((psiti(valid)-potential)/1000+total*atr)./resistance(valid);
    [lowest,index]=min(atrani);
    if lowest>=-1e-12, atrani=max(0,atrani); atr=sum(atrani); return; end
    flag(index)=true;
end
error('OFM:Roots','Root uptake allocation failed.');
end

function q=vert(k1,k2,ksat1,ksat2,thick1,thick2,psi1,psi2,stone1,stone2,rhowg)
if k1<=0 || k2<=0, q=0; return; end
meanK=exp((thick2*log(k1)+thick1*log(k2))/(thick1+thick2));
meanK=min([meanK ksat1 ksat2]);
q=(psi1-psi2)/((thick1+thick2)/2)*meanK/rhowg*(1-(stone1+stone2)/2);
end
