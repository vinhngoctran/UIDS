%============================================================================
%   soilground.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program simulates daily transpiration and soil water fluxes through 
%   a soil profile covered with vegetation. This model is developed based 
%   on the concept of Soil Vegetation Atmosphere Transport (SVAT) model 
%   (Hammel & Kennel, 2001) written in Fortran.
%   (https://www.lwf.bayern.de/boden-klima/wasserhaushalt/index.php)
%============================================================================

function [psim,gwat,vrfli,infli,byfli,dsfli,ntfli,trani,swat,gwfl,seep] = soilground...
                                            (modsfl_inf,psim,ptr,gwat,slvp,slfl,deltaT)

%% Parameters and characteristics
nlayer = modsfl_inf.nlayer;             
thick = modsfl_inf.thick;                        
stonef = modsfl_inf.stonef;         
psicr = modsfl_inf.psicr;            
par = modsfl_inf.par;           
ilayer = modsfl_inf.ilayer;           
infexp = modsfl_inf.infexp;           
qlayer = modsfl_inf.qlayer;            
dispc = modsfl_inf.dispc;             
height = modsfl_inf.height;                                       
lai = modsfl_inf.lai;                                     
sai = modsfl_inf.sai;             
snow = modsfl_inf.snow;              
snoden = modsfl_inf.snoden;           
mxrtln = modsfl_inf.mxrtln;          
mxkpl = modsfl_inf.mxkpl;          
densef = modsfl_inf.densef;         
frelden = modsfl_inf.frelden;            
tini = modsfl_inf.tini;              
age = modsfl_inf.age;               
rgroper = modsfl_inf.rgroper;          
inirdep = modsfl_inf.inirdep;         
inirlen = modsfl_inf.inirlen;         
rtrad = modsfl_inf.rtrad;             
fxylem = modsfl_inf.fxylem;          
qffc = modsfl_inf.qffc;            
qfpar = modsfl_inf.qfpar;         
dslope = modsfl_inf.dslope;         
length = modsfl_inf.length;           
drain = modsfl_inf.drain;          
gsc = modsfl_inf.gsc;            
gsp = modsfl_inf.gsp;             

% Example
% nlayer = 1;             % The number of soil layers
% thick = 2;              % Thickness of each soil layers [mm]
% stonef = 0.5;           % stone volume fraction, unitless
% psicr = 50;             % minimum plant leaf water potential, MPa
% par(1,1) = 5;           % theta at saturation,                                   
% par(10,1) = 0.3;        % wetness at dry end of near-saturation range
% par(7,1) = 0.01;        % clapp and hornberger m, kpa
% par(8,1) = 0.001;       % clapp and hornberger m, kpa
% par(6,1) = 0.01;        % saturated hydraulic conductivity
% par(9,1) = 0.1;         % exponent for psi-theta relation
% par(2,1) = 5;           % volumetric water content at field capacity
% par(4,1) = 0.1;         % matric potential at field capacity, kpa
% par(5,1) = 0.1;         % wetness at field capacity, dimensionless
% par(3,1) = 2;           % hydraulic conductivity at field capacity, mm/d
% ilayer = 1;             % number of layers over which infiltration is distributed
% infexp = 0;             % infiltration exponent, 0 all to top, 1 uniform with depth, >1.0=more at bottom than at top
% qlayer = 0;             % number of soil layers for SRFL           
% dispc = 1;              % zero-plane displacement for closed canopy, m
% height = 2;             % canopy height, m, minimum of 0.01 m                                  
% lai = 0.2;              % leaf area index, m2/m2, minimum of 0.00001                           
% sai = 0.1;              % stem area index, m2/m2
% snow = 0;               % water equivalent of snow on the ground, mm
% snoden = 0;             % snow density, mm/mm
% mxrtln = 0.1;           % maximum root length per unit land area, m/m2
% mxkpl = 0.01;           % maximum plant conductivity, (mm/d)/MPa
% densef = 0.01;          % density factor
% frelden = 1;            % final relative values of root length per unit volume
% tini = 0;               % initial time for root growth in layer
% age = 1;                % age of vegetation
% rgroper = 1;            % period of root growth in layer, a
% inirdep = 0.5;          % intial root depth, m
% inirlen = 0.3;          % intial total root length, m m-2
% rtrad = 10;             % average root radius, mm
% fxylem = 0.1;           % fraction of plant resistance in xylem
% qffc = 0.5;             % BYFL fraction at field capacity
% qfpar = 0.001;          % quick flow parameter
% dslope = 2;             % slope for soil water flow, radians; no DSFLI if DSLOPE = 0
% length = 1;             % slope length (for DSFLI), m
% drain = 0.01;           % gravity drainage (should very small)
% gsc = 0.2;              % discharge from GWAT, fraction per day, d-1
% gsp = 0.2;              % fraction of discharge to seepage        
            
%% States
% Convert from m/s to mm/day
dt = deltaT/86400;                          % Computational timestep
ptr = ptr * 1000 * 86400;                   % average potential transpiration rate over time period, mm/d       RESULTS FROM CANOPY/EB MODEL
slvp = slvp * 1000 * 86400;                 % evaporation rate from soil, mm/d                                  RESULS FROM CANOPY OR EB MODEL
slfl = slfl * 1000 * 86400;                 % input rate to soil surface, mm/d                                  RESULTS FROM INFILTRATION
gwat = gwat;                                % groundwater storage below soil layers, mm                         INITIAL CONDITION and UPDATE CONTINUOUSLY
psim = psim;                                % matric soil water potential, kPa                                  INITIAL STATE and UPDATE CONTINUOUSLY
rhowg = 0.0098;                             % density of water times gravity acceleration, MPa/m or kPa/mm

% Example
% dt = 1;
% ptr = 10; 
% slvp = 1;    
% slfl = 10; 
% gwat = 0.1; 
% psim = 500;
%% MODEL SIMULATION
[psig, swatmx, wetf, wetc, chm, chn, wetnes, swati, ksat, par] =  soilpar(nlayer, par, thick, stonef, psim, psicr);
[psiti, theta, swat, kk] = soilvar(nlayer, par, psig, psim, wetnes, swati);
infrac = infpar(infexp, ilayer, nlayer, thick);
% [swatqx, swatqf] = srfpar(qlayer, par, thick, stonef, swatmx);
relden = rootgrowth (frelden, tini, age, rgroper, inirdep, inirlen, nlayer);
[height, lai, sai, rtlen, rplant] =  canopy (snow, snoden, mxrtln, mxkpl, densef, height, lai, sai);
[rxylem, rrooti, alpha] = plntres (nlayer, thick, stonef, rtlen, relden, rtrad, rplant, fxylem);
[atr,atrani] = tbylayer(ptr, dispc, alpha, kk, rrooti, rxylem,psiti, nlayer, psicr);
for i = 1:nlayer
    atri(i) = atrani(i);
end
% if (atr < ptr)
%     ger = swge(aa, asubs, vpd, raa, ras, rss, delta, atr);
% end
for i = 1: nlayer
     trani(i) = atri(i) / dt;
end
byfrac = byflfr(nlayer, wetnes, par, qffc, qfpar);
for i = 1:nlayer
    k = nlayer+1-i;
    % downslope flow rates
    dsfli = dslop(dslope, length, thick, stonef, psim, rhowg, kk);
    % vertical flow rates
    if (i < nlayer)
        if (abs(psiti(i) - psiti(i + 1)) < dpsimx)
            vrfli(i) = 0.0d0;
        else
            vrfli = vert(kk, kk1, ksat, ksat1, thick, thick1, psit, psit1, stone, stone1, rhowg);
        end
    else
        %  bottom layer
        if (drain > 0.00001)
            %  gravity drainage only
            vrfli(nlayer) = drain * kk(nlayer) * (1.0d0 - stonef(nlayer));
        else
            vrfli(nlayer) = 0.0d0;
        end
    end
end
% net inflow to each layer including E and T withdrawal adjusted for interception
[vv, byfli, infli, ntfli] = inflow(nlayer, dt, infrac, byfrac, slfl, vrfli, dsfli, trani, slvp, swatmx,swati);
for i = 1: nlayer
    vrfli(i) = vv(i);
end
% groundwater flow and seepage loss
[gwfl, seep] = gwater(gwat, gsc, gsp, dt, vrfli(nlayer));
%  end of rate calculations
%  integrate below ground storages over iteration interval
for i = 1:nlayer
    swati(i) = swati(i) + ntfli(i) * dt;
end
gwat = gwat + (vrfli(nlayer) - gwfl - seep) * dt;         %! groundwater storage below soil layers, mm
% new soil water variables and test for errors
for i = 1: nlayer
    ths=par(1,i);
    thr=par(10,i);
    wetnes(i) = (ths * swati(i) / swatmx(i) -thr) / (ths - thr);
    if wetnes(i) < 0
        wetnes(i) = 0;
    end
    psim(i) = fpsim(wetnes(i),par(:,i));
end
[psiti, theta, swat, kk] = soilvar (nlayer, par, psig, psim, wetnes, swati);
% display('OK check')
% Check balance error;
BalanceE = (slfl - swat/dt - gwat/dt)/slfl*100;

%% Main Results
% vrfli           %! vertical drainage rate from layer i, mm/d
% infli           %! infiltration rate into layer, mm/d
% byfli           %! bypass flow rate from layer, mm/d
% dsfli           %! downslope flow rate from layer, m/d
% ntfli           %! net flow rate into layer, mm/d
% trani           %! transpiration rate from layer, mm/d
% swat            %! total soil water in all layers, mm
% gwfl            %! streamflow from groundwater discharge, mm/d
% seep            %! deep seepage loss from groundwater, mm/d
% gwat            %! groundwater storage below soil layers, mm

end

%%
%============================================================================
%============================================================================
%============================================================================
%   soilground.m
%
%   Project:    PCUW
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Sub-Program simulates daily transpiration and soil water fluxes through 
%   a soil profile covered with vegetation. This model is developed based 
%   on the concept of Soil Vegetation Atmosphere Transport (SVAT) model 
%   (Hammel & Kennel, 2001) written in Fortran.
%   (https://www.lwf.bayern.de/boden-klima/wasserhaushalt/index.php)
%============================================================================
%============================================================================
%============================================================================
%============================================================================
%%
function [psig, swatmx, wetf, wetc, chm, chn, wetnes, swati, ksat, par] =  soilpar...
    (nlayer, par, thick, stonef, psim, psicr)

% calculated soil water parameters and initial variables
% Input variables
% nlayer     number of soil layers
% thick(*)   layer thicknesses, mm"
% thsat      theta at saturation, matrix porosity, imodel=0
% ths        theta at saturation, matrix porosity, imodel=1
% thr        residual theta, imodel=1
% stonef(*)  stone volume fraction, unitless"
% thetaf     volumetric water content at field capacity"
% psif       matric potential at field capacity, kpa
% bexp       exponent for psi-theta relation
% wetinf     wetness at dry end of near-saturation range
% psim(*)    matric soil water potential for layer, kpa
% kf         hydraulic conductivity at field capacity, mm/d
% psicr      minimum plant leaf water potential, mpa

% output
% psig(*)    gravity potential, kpa
% swatmx(*)  maximum water storage for layer, mm
% wetf       wetness at field capacity, dimensionless
% wetc(*)    wetness at psicr, dimensionless
% chm        clapp and hornberger m, kpa
% chn        clapp and hornberger n
% wetnes(*)  wetness, fraction of saturation
% swati(*)   water volume in layer, mm
% ksat       saturated hydraulic conductivity, mm/d
rhowg = 0.00981d0 ;   % RHOWG  density of water times gravity acceleration, MPa/m or kPa/mm
for i = 1: nlayer
    % gravity potential is negative down from surface
    if (i == 1)
        psig(1) = -rhowg * thick(1) / 2.0d0;
    else
        psig(i) = psig(i - 1) - rhowg * ((thick(i - 1) + thick(i)) / 2.0d0);
    end
    thsat = par(1,i);
    thetaf=par(2,i);
    swatmx(i) = thick(i) * thsat * (1.0d0 - stonef(i));
    wetf = thetaf / thsat;
    par(5,i) = wetf;
    wetinf=par(10,i);
    bexp=par(9,i);
    psif=par(4,i);
    psiinf(i) = psif * (wetinf / wetf)^(-bexp);
    chm = (-psiinf(i) / (1.0d0 - wetinf)^2) - bexp * (-psiinf(i)) / (wetinf * (1.0d0 - wetinf));
    par(7,i) = chm;
    chn = 2.0d0 * wetinf - 1.0d0 - (-psiinf(i) * bexp / (chm * wetinf));
    par(8,i) = chn;
    
    if (psim(i) >= 0.0d0)           %print*, 'matrix psi must be negative or zero'
        wetnes(i) = 1.0d0;
    else
        wetnes(i) = wetf * (psim(i) / psif)^(-1.0d0 / bexp);
        if (wetnes(i) > wetinf)
            wetnes(i) = (1.0d0 + chn) / 2.0d0 + 0.50d0 * sqrt(chn^2 - 2.0d0 * chn + 1.0d0 + 4.0d0 * psim(i) / chm);
        end
    end
    swati(i) = wetnes(i) * swatmx(i);
    kf=par(3,i);
    ksat = kf * (1.0d0 / wetf)^(2.0d0 * bexp + 3);
    par(6,i) = ksat;
    wetc(i) = wetf * (1000 * psicr / psif)^(-1.0d0 / bexp);
end
end

%%
function [psiti, theta, swat, kk] = soilvar (nlayer, par, psig, psim, wetnes, swati)

% soil water variables

% Input variables
% nlayer     number of soil layers
% mpar,ml    maximum number of parameters and layers
% psig(*)    gravity potential, kpa
% psim(*)    matric soil water potential for layer, kpa
% wetnes(*)  wetness, fraction of saturation.
% thsat      theta at saturation, matrix porosity
% kf         hydraulic conductivity at field capacity, mm/d
% bexp       exponent for psi-theta relation
% wetf       wetness at field capacity, dimensionless
% swati(*)   water volume in layer, mm

% output
% psiti(*)   total potential, kpa
% theta(*)   water content, mm water / mm soil matrix
% swat       total soil water in all layers, mm
% kk(*)      hydraulic conductivity, mm/d

swat = 0.0d0;
for i = 1: nlayer
    kf = par(3,i);
    wetf = par(5,i);
    bexp=par(9,i);
    psiti(i) = psim(i) + psig(i);
    theta(i)=wetnes(i)*par(1,i);  
    if wetnes(i) > 0.0001
        kk(i)= kf * (wetnes(i)/wetf).^(2*bexp+ 3.0d0);
    else
        %            extremely dry
        kk(i) = 1e-10;
    end    
    swat = swat + swati(i);
end
end

%%
function [vv, byfli, infli, ntfli] = inflow(nlayer, dti, infrac, byfrac, slfl, vrfli, dsfli, trani, slvp, swatmx,swati)
% Inflow and byflow for each layer, and net inflow including e and t withdrawal

% Input variables
%     nlayer     number of soil layers being used, max of 20
%     dti        time step for iteration interval, d
%     infrac(*)  fraction of infiltration to each layer
%     byfrac(*)  fraction of layer infiltration to bypass flow
%     slfl       input rate to soil surface, mm/d
%     dsfli(*)   downslope flow rate from layer, mm/d
%     trani(*)   transpiration rate from layer, mm/d
%     slvp       evaporation rate from soil, mm/d
%     swatmx(*)  maximum water storage for layer, mm
%     swati(*)   water volume in layer, mm
%     vrfli(*)   vertical drainage rate from layer, mm/d

% output
%     vv(*)      modified vrfli, mm/d
%     byfli(*)   bypass flow rate from layer, mm/d
%     infli(*)   infiltration rate into layer, mm/d
%     ntfli(*)   net flow rate into layer, mm/d

    for k = 1: nlayer
        i = nlayer+1-k;
%        need to do from bottom up
        infil = slfl * infrac(i);
%         pause
        byfli(i) = byfrac(i) * infil;
%         pause
        infli(i) = infil - byfli(i);
%         pause
        if (i == nlayer)
            vv(i) = vrfli(i);
        end
        if (i > 1)
            maxin = (swatmx(i) - swati(i)) / dti + vv(i) + dsfli(i) + trani(i);
            if (vrfli(i - 1) + infli(i) > maxin)
%               adjust to prevent oversaturation
                if (byfrac(i) > 0)
                    if (vrfli(i - 1) > maxin) 
%                    reduce infli, increase byfli
                        byfli(i) = byfli(i) + infli(i) - (maxin - vrfli(i - 1));
                        infli(i) = maxin - vrfli(i - 1);
                        vv(i - 1) = vrfli(i - 1);
                    else
%                    shift all infli to byfli and reduce vrfli(i - 1)
                        byfli(i) = byfli(i) + infli(i);
                        infli(i) = 0;
                        vv(i - 1) = maxin;
                    end
                else
%                 no bypass flow allowed, reduce vrfli(i-1), infli(i) unchanged
                    vv(i - 1) = maxin - infli(i);
                end
            else
%               byfli and infli unchanged
                vv(i - 1) = vrfli(i - 1);
            end
            ntfli(i) = vv(i - 1) + infli(i) - vv(i) - dsfli(i) - trani(i);
        else
%            i = 1
            maxin = (swatmx(1) - swati(1)) / dti + vv(1) + dsfli(1) + trani(1) + slvp;
%             pause
            if (infli(1) > maxin) 
%               increase byfli(1) to prevent oversaturation
                byfli(1) = byfli(1) + infli(1) - maxin;
                infli(1) = maxin;
%                  may be negative
            end
            ntfli(1) = infli(1) - vv(1) - dsfli(1) - trani(1) - slvp;
        end
    end
end


%%
function  byfrac = byflfr(nlayer, wetnes, par, qffc, qfpar)

% bypass flow fraction of infiltration to layer

% input
%     bypar      1 to allow byfl, or 0 to prevent byfl (bypass flow)
%     nlayer     number of soil layers to be used in model, <= ml
%     wetnes(*)  wetness, fraction of saturation
%     wetf(*)    wetness at field capacity, dimensionless
%     qffc       byfl fraction at field capacity
%     qfpar      quick flow parameter
    
% output
%     byfrac(*) ! fraction of layer infiltration to bypass flow

bypar = 1; 
    for i = 1: nlayer
        if (bypar >= 1) 
            byfrac(i) = qffc^(1.0d0 - (1.0d0 / qfpar) * (wetnes(i) - par(5,i)) / (1.0d0 - par(5,i)));
            if (byfrac(i) > 1)
                byfrac(i) = 1;
            end
%            generate bypass flow to avoid saturation
            if (wetnes(i) > 0.99)
                byfrac(i) = 1;
            end
        else
            byfrac(i) = 0;
        end
    end
    
end

%%
function [height, lai, sai, rtlen, rplant] =  canopy (snow, snoden, mxrtln, mxkpl, densef, height, lai, sai)
%  canopy parameters

% input
%     snow       water equivalent of snow on the ground, mm
%     snoden     snow density, mm/mm
%     mxrtln     maximum root length per unit land area, m/m2
%     mxkpl      maximum plant conductivity, (mm/d)/mpa
%     cs         ratio of projected sai to canopy height, m-1, not needed in this version
%     densef     density factor

% output
%     height     canopy height above any snow, m, minimum of 0.01 m
%     lai        leaf area index, m2/m2, minimum of 0.00001
%     sai        stem area index, m2/m2
%     rtlen      root length per unit land area, m/m2
%     rplant     plant resistivity to water flow, mpa d/mm


    snodep = .0010d0 * snow / snoden;
    hnosno = max(.010d0, height);
    hsno = max(0.0d0, hnosno - snodep);
    ratio = hsno / hnosno;
    height = max(.010d0, hsno);

    lai = ratio * densef * lai;
%       sai = densef * cs * height           original brook90 expression replaced by direct input via climate.in
    sai = densef * sai;
    if (lai < .000010d0) 
        lai = .000010d0;
    end

    rtlen = densef * mxrtln;
    kpl = densef * mxkpl;
    if (kpl < 1e-08) 
        kpl = 1e-08;
    end
    rplant = 1.0d0 / kpl;

end

%%
function dsfli = dslop(dslope, length, thick, stonef, psim, rhowg, kk)
% downslope flow rate from layer

% input
%     dslope   slope for soil water flow, radians
%                          ! no dsfli if dslope = 0
%     length   slope length (for dsfli), m
%     thick    layer thicknesses, mm
%     stonef   stone volume fraction, unitless
%     psim     matric soil water potential, kpa
%     rhowg    density of water times acceleration of gravity,kpa/mm
%     kk       hydraulic conductivity, mm/d

% output
%     dsfli    downslope flow rate from layer, mm/d
 
    ll = 1000 * length;
    grad = rhowg * sin(dslope) + (2 * psim / ll) * cos(dslope);
%       grad = rhowg * sin(dslope)
    aratio = thick * (1.0d0 - stonef) * cos(dslope) / ll;
    dsfli = kk * aratio * grad / rhowg;
%     no water uptake into dry soil because no free water at outflow face
    if dsfli < 0
        dsfli = 0.0d0;
    end
end

%%
function fk = fk(wetnes, par, imodel)
% hydraulic conductivity from wetness

% input
%     wetnes    wetness, fraction of saturation
%     imodel    parameterization of hydraulic functions
%     psif      matrix potential at field capacity, kpa
%     bexp      exponent for psi-theta relation
%     wetinf    wetness at dry end of near-saturation range
%     wetf      saturation fraction at field capacity
%     chm       clapp and hornberger m, kpa
%     chn       clapp and hornberger n
%     kf        hydraulic conductivity at field capacity, mm/d
%     ths       water content at saturation
%     thr       residual water content
%     alfa      parameter alpha, 1/m
%     mvgn      parameter n
%     ks        hydraulic conductivity at saturation, mm/d
%     a         tortuosity parameter (default = 0.5)

% output
%     fk         hydraulic conductivity, mm/d

    fk = 0.0d0;
 
    if (imodel == 0)
        kf=par(3);
        wetf=par(5);
        bexp=par(9);
        if (wetnes > 0.0001) 
            fk=kf * (wetnes / wetf)^(2.0d0 * bexp+ 3.0d0);
        else
% !      extremely dry
            fk = 1e-10;
        end
    end
    if (imodel == 1) 
        ks=par(6);
        alfa=par(7);
        mvgn=par(8);
        a=par(9);
        awet=max(wetnes,1.e-6);
        awet=min(wetnes,1.0d0);
        fk =ks*awet^a*(1.0d0-(1.0d0-awet^(mvgn/(mvgn-1.0d0)))*(1.0d0-1.0d0/mvgn))^2;
% !        write(10,*) fk,awet,log(fk)
    end
end

%%
function fpsim = fpsim(wetnes, par)
% matric potential from wetness

% input
%     wetnes    wetness, fraction of saturation
%     imodel    parameterization of hydraulic functions
%     par(*)    parameter array
%     psif      matrix potential at field capacity, kpa
%     bexp      exponent for psi-theta relation
%     wetinf    wetness at dry end of near-saturation range
%     wetf      saturation fraction at field capacity
%     chm       clapp and hornberger m, kpa
%     chn       clapp and hornberger n
%     ths       water content at saturation
% 	  thr       residual water content
%     alfa      parameter alpha, 1/m
%     mvgn      parameter n
%     mvgm      parameter m
%     ks        hydraulic conductivity at saturation, mm/d
%     a         tortuosity parameter (default = 0.5)

% output
%     fpsim     matric potential, kpa

 imodel = 1;
    fpsim = 0.0d0;
 
    if(imodel == 0)
        psif=par(4);
        bexp=par(9);
        wetf=par(5);
        wetinf=par(10);
        chm=par(7);
        chn=par(8);
        if (wetnes <= 0.0d0) 
%        arbitrary very negative value
            fpsim = -1e+10;
        elseif (wetnes < wetinf)
            fpsim = psif * (wetnes / wetf)^(-bexp);
        elseif (wetnes < 1.0d0)
%         in near-saturated range
            fpsim = chm * (wetnes - chn) * (wetnes - 1.0d0);
        else
%         saturated
            fpsim = 0.0d0;
        end
    end
    if (imodel == 1)
        alfa=par(7);
        mvgn=par(8);
        mvgm=1.0d0-1.0d0/mvgn;
        awet=max(wetnes,1.e-6);
        awet=min(wetnes,1.0d0);
        fpsim = (-1.0d0/alfa)*(awet^(-1.0d0/mvgm)-1)^(1.0d0/mvgn);
        fpsim = fpsim * 9.810d0;  %! conversion from m to kpa
    end
end

%%
function [gwfl, seep] = gwater(gwat, gsc, gsp, dt, vrfln)
% calculates groundwater flow and seepage loss

% input
%     gwat     groundwater storage below soil layers, mm
%     gsc      discharge from gwat, fraction per day, d-1
%     gsp      fraction of discharge to seepage
%     dt       time step for interval, d
%     vrfln    vertical drainage rate from lowest layer, mm/d

% output
%     gwfl     streamflow from groundwater discharge, mm/d
%     seep     deep seepage loss (baseflow) from groundwater, mm/d
 
    if (gsc < 1e-08)
%         no groundwater
        seep = gsp * vrfln;
        gwfl = vrfln - seep;
    else
        seep = gwat * gsc * gsp;
        gwfl = gwat * gsc * (1. - gsp);
%         prevent negative gwat
        if (gwat / dt - (gwfl + seep) < 0)
            seep = gsp * gwat / dt;
            gwfl = (1 - gsp) * gwat / dt;
        end
    end
end

%%
function  infrac = infpar(infexp, ilayer, nlayer, thick)

% input
%     infexp     infiltration exponent, 0 all to top, 1 uniform with depth
%                                 >1.0=more at bottom than at top
%     ilayer     number of layers over which infiltration is distributed
%     nlayer     number of soil layers being used
%     thick(*)   layer thicknesses, mm
% output
%     infrac(*)  fraction of infiltration to each layer
 
    if (infexp <= 0.)
        infrac(1) = 1.0d0;
        for i = 2:nlayer
            infrac(i) = 0.0d0;
        end
    else
        thickt = 0.0d0;
        for i = 1: ilayer
            thickt = thickt + thick(i);
        end
        thicka = 0.0d0;
        for i = 1: nlayer
            if (i <= ilayer)
                if i==1
                    pthick = thicka;
                    thicka(i) = pthick + thick(i);
                    infrac(i) = (thicka(i) / thickt)^infexp - (pthick / thickt)^infexp;
                else
                    thicka(i) = thicka(i - 1) + thick(i);
                    infrac(i) = (thicka(i) / thickt)^infexp - (thicka(i - 1) / thickt)^infexp;
                end
            else
                infrac(i) = 0.0d0;
            end
        end
    end
end 

%%
function interp = interp (npairs, funct, xvalue,xx,yy)
%      interpolates between points in data functions

% nput
%     npairs   number of pairs of values to be used
%     funct(*) array of pairs of values: x1, y1, x2, y2, ...
%     xvalue   x value
% output
%     interp   y value


interp = 0;
%      put funct into xx and yy
i = 0;
id = [1:2:2 * npairs - 1];
for j = 1: numel(id)
    
    i = i + 1;
    xx(i) = funct(id(j));
    yy(i) = funct(id(j) + 1);
end
%      interpolate using xx and yy
for j = 1: npairs
    if (xvalue == xx(j))
        interp = yy(j)
        
    elseif (xvalue < xx(j))
        interp = yy(j - 1) + (xvalue - xx(j - 1)) * (yy(j) - yy(j - 1)) / (xx(j) - xx(j - 1))
  
    end
end
 
end
 
%%
function dtinew = iter(nlayer, dti, dpsidw, ntfli, swatmx, psiti, dswmax, dpsimx, thick, wetnes, par, imodel, trani,slvp, mpar, ml, dtimin, pr)

% input
%     nlayer     number of soil layers to be used in model
%     dti        time step for iteration interval, d
%     dtimin     minimum time step for iteration interval, d
%     dpsidw(*)  rate of change of total potential with water content, kpa/mm
%     ntfli(*)   net flow rate into layer, mm/d
%     thick(*)   thickness of layer, mm
%     wetnes(*)  water saturation of layer
%     trani(*)   actual transpiration (also output after check)
%     slvp       actual evaporation (also output after check)
%     swatmx(*)  maximum water storage for layer, mm
%     swati(*)   actual water storage for layer, mm
%     psiti(*)   total potential, kpa
%     dswmax     maximum change allowed in swati, percent of swatmx(i)
%     dpsimx     maximum potential difference considered "equal", kpa  

% output
%     dtinew     second estimate of dti

%      first approximation to new total potential
    for i = 1: nlayer
        if (imodel == 0)
            a(i) = ntfli(i) * dpsidw(i) / swatmx(i);
        end
        if (imodel == 1)
            a(i) = ntfli(i)/thick(i) * dpsidw(i) / (par(1,i) - par(10,i));
        end
        t(i) = psiti(i) + a(i) * dti;
    end
%      test to see if dti should be reduced
    dtinew = dti;
    for i = 1:nlayer
%         prevent too large a change in water content
        dtinew = min(dtinew, .010d0 * dswmax * swatmx(i) / max(.0000010d0, abs(ntfli(i))));
%         prevent a change in water content larger than total available water
        if((imodel == 1) && (ntfli(i) < 0))
            thr=par(10,i);
            th=wetnes(i)*par(1,i);
            dtinew = min(dtinew,(thr-th)*thick(i)/ntfli(i)/1.30d0);
            if(dtinew < dtimin) 
                thr=par(10,i);
                
                th=wetnes(i)*par(1,i);
                psi=fpsim(wetnes(i),par(1,i),imodel);
                k= fk(wetnes(i),par(1,i),imodel);
                for j = 1: nlayer
                    th=wetnes(j)*par(1,j);
                    thr=par(10,j);
                    psi=fpsim(wetnes(j),par(1,j),imodel);
                    k= fk(wetnes(j),par(1,j),imodel);
 
                    if (pr == 1) 
%                         call realpr('xxx i=, th=, thr=, netflow=, thick=, k=, psi=', -1,(/real(j,8),th,thr,ntfli(j),thick(j),k,psi/), 7)
                    end
                end
                dtinew=dtimin;
                trani(i)=0;
                if(i == 1) 
                    slvp=0;
                end
            end
        end
        if((imodel == 0) && (ntfli(i)< 0))
            th=wetnes(i)*par(1,i);
            dtinew = min(dtinew,-th*thick(i)/ntfli(i)/1.30d0);
            if (dtinew < dtimin) 
                for j = 1: nlayer
                    th=wetnes(j)*par(1,j);
                    thr=par(10,j);
                    psi=fpsim(wetnes(j),par(1,j),imodel);
                    k= fk(wetnes(j),par(1,j),imodel);
 
                    if (pr == 1)
%                         call realpr('xxx i=, th=, netflow=, thick=, k=, psi=', -1, (/real(j,8),th,ntfli(j),thick(j),k,psi/),7)
                    end
                end
                dtinew=dtimin;
                trani(i)=0;
                if (i == 1) 
                    slvp=0;
                end
            end
        end
%         prevent oscillation of potential gradient
        if (imodel == 0) 
            if (i < nlayer) 
%            total potential difference at beginning of iteration
                pp = psiti(i) - psiti(i + 1);
%            first approx to total potential difference at end of iteration
                tt = t(i) - t(i + 1);
                if (abs(tt) > dpsimx && abs(pp) < dpsimx) 
                    if (tt < 0.0d0 && pp > 0.0d0 || tt > 0.0d0 && pp < 0.0d0) 
                        dtinew = min(dtinew, -pp / (a(i) - a(i + 1)));
                        dtinew = max(dtinew, dtimin);
                    end
                end
            end
        end
    end
end

%%

function [rxylem, rrooti, alpha] = plntres (nlayer, thick, stonef, rtlen, relden, rtrad, rplant, fxylem)
% allocates total plant resistance to xylem and root layers

% input
%     nlayer    number of soil layers (max 50)
%     thick(*)  layer thicknesses, mm
%     stonef(*) stone volume fraction, unitless
%     rtlen     root length per unit land area, m/m2,
%     relden(*) relative values of root length per unit volume
%     rtrad     average root radius, mm
%     rplant    plant resistance to water flow, mpa d/mm,
%     fxylem    fraction of plant resistance in xylem

% output
%     rxylem    xylem resistance, mpa d/mm, 1e20 if no roots
%     rrooti(*) root resistance for layer, mpa d/mm, 1e20 if no roots
%     alpha(*)  modified cowan alpha, mpa


rhowg = 0.00981d0 ;   %! RHOWG  density of water times gravity acceleration, MPa/m or kPa/mm
    rxylem = fxylem * rplant;

    suM = 0.0d0;
    for i = 1:nlayer
        d(i) = thick(i) * (1.0d0 - stonef(i));
        suM = suM + relden(i) * d(i);
    end
    for i = 1:nlayer
        if (relden(i) < 0.000010d0 || rtlen < 0.10d0) 
%            no roots in layer
            rrooti(i) = 1e+20;
            alpha(i) = 1e+20;
        else
            rtfrac = relden(i) * d(i) / suM;
%            root resistance for layer
            rrooti(i) = (rplant - rxylem) / rtfrac;
%            rhizosphere resistance for layer
            rtdeni = rtfrac * .0010d0 * rtlen / d(i);
%                             .001 is (mm/mm2)/(m/m2) conversion
            delt = pi * rtrad^2.0d0 * rtdeni;
            alpha(i) = (1.0d0 / (8.0d0 * pi * rtdeni)) * (delt - 3.0d0 - 2.0d0 * (log(delt)) / (1.0d0 - delt));
            alpha(i) = alpha(i) * .0010d0 * rhowg / d(i);
%                             .001 is mpa/kpa conversion
        end
    end

end


%%
function  relden = rootgrowth (frelden, tini, age, rgroper, inirdep, inirlen, nlayer)
% root growth according to lwf root growth model, hammel and kennel, 2000

% input
%     frelden(*) final relative values of root length per unit volume
%     tini(*)    initial time for root growth in layer
%     age        age of vegetation
%     rgroper    period of root growth in layer, a
%     inirdep    intial root depth, m
%     inirlen    intial total root length, m m-2
%     rl0        constant intial root length density, m m-3
%     nlayer     number of soil layers (max 1000)
% 
% output        
%    relden(*)   actual relative values of root length per unit volume
    if(rgroper > 0.0d0) 
        for i = 1: nlayer
            if(age < tini(i)) 
                relden(i)=0.0d0;
            end
            if ((age >= tini(i)) && (age <= rgroper+tini(i))) 
                rl0=inirlen/inirdep;
                relden(i)=rl0*(frelden(i)/rl0)^((age-tini(i))/rgroper);
            end 
            if (age < rgroper+tini(i)) 
                relden(i)=frelden(i);
            end
        end
    else
        for i = 1:nlayer
            relden(i)=frelden(i);
        end
    end

end

%%
function [swatqx swatqf] = srfpar (qlayer, par, thick, stonef, swatmx)
% source area parameters
% input
%     qlayer     number of soil layers for srfl
%     mpar,ml    maximum number of parameters and layers
%     thetaf(*)  volumetric water content of layer at field capacity
%     thick(*)   layer thickness, mm
%     stonef(*)  stone volume fraction of layer
%     swatmx(*)  maximum water storage for layer, mm

% output
%     swatqx     maximum water storage for layers 1 through qlayer, mm
%     swatqf     water storage at field capacity for layers 1 through qlayer, mm

    swatqx = 0.0d0;
    swatqf = 0.0d0;
    for i = 1: qlayer
        swatqx = swatqx + swatmx(i);
        swatqf = swatqf + par(2,i) * thick(i) * (1.0d0 - stonef(i));
    end
 
end

%%
function  erate = swge(aa, asubs, vpd, raa, ras, rss, delta, arate)
%     shuttleworth and wallace (1985) ground evaporation when transpiration known
% input
%     aa       net radiation at canopy top minus ground flux, w/m2
%     asubs    net radiation minus ground flux at ground, w/m2
%     vpd      vapor pressure deficit, kpa
%     raa      boundary layer resistance, s/m
%     ras      ground-air resitance, s/m
%     rss      ground evaporation resistance, s/m
%     delta    desat/dtair, kpa/k
%     arate    actual transpiration rate, mm/d
% output
%     erate    ground evaporation rate, mm/d

    lec = arate / (etom * wtomj);
    rs = (delta + gamma) * ras + gamma * rss;
    ra = (delta + gamma) * raa;
    le = (rs / (rs + ra)) * lec + (cprho * vpd + delta * ras * asubs + delta * raa * aa) / (rs + ra);
    erate = etom * wtomj * (le - lec);

end

%%
function [atr,atrani] = tbylayer(ptr, dispc, alpha, kk, rrooti, rxylem,psiti, nlayer, psicr)

% actual transpiration rate by layers

% Input variables
%       
% j         1 for daytime, 2 for nighttime
% ptr       average potential transpiration rate over time period, mm/d
% dispc     zero-plane displacement for closed canopy, m
% alpha     modified cowan alpha, mpa
% kk	    hydraulic conductivity, mm/d
% rrooti    root resistance for layer, mpa d/mm
% rxylem    xylem resistance, mpa d/mm
% psiti     total soil water potential, kpa
% nlayer    number of soil layers (max 20)
% psicr     critical potential for plant, mpa
% nooutf   	1 if no outflow allowed from roots, otherwise 0

% output
% atr       actual transpiration rate over time period, mm/d
% atrani	actual transpiration rate from layer over time period, mm/d

rhowg = .00981d0; %    RHOWG  density of water times gravity acceleration, MPa/m or kPa/mm
atr = 0;
atrani = 0;
j = 1;
nooutf = 1;
for i = 1: nlayer
    if (rrooti(i) > 1e+15)
        flag(i) = 1;
    end
    if (nooutf == 1 && psiti(i) / 1000.0d0 <= psicr)
        flag(i) = 1;
    else
        flag(i) = 0;
    end
end

% top of loop for recalculation of transpiration if more layers get flagged
for ii = 1:100
    negflag = 0;
    suM = 0.0d0;
    for i = 1:nlayer
        if (flag(i) == 0)
            ri(i) = rrooti(i) + alpha(i) / kk(i);
            suM = suM + 1.0d0 / ri(i);
        else
            %               flagged
            atrani(i) = 0.0d0;
        end
    end
    if (suM < 1e-20)
        %           all layers flagged, no transpiration
        atr = 0.0d0;
        psit = -1e+10;
    else
        rt = 1.0d0 / suM;
        
        %        weighted mean soil water potential
        psit = 0.0d0;
        for i = 1:nlayer
            if (flag(i) == 0)
                psit = psit + rt * psiti(i) / ri(i);
            end
        end
        %       soil water supply rate, assumed constant over day
        supply = (psit / 1000.0d0 - psicr - rhowg * dispc) / (rt + rxylem);
        %        transpiration rate limited by either ptr or supply
        if (j == 1)
            %            daytime, ptr is average of a half sine over daytime
            r = (2.0d0 / pi) * (supply / ptr);
            if (r <= 0.0d0)
                atr = 0.0d0;
            elseif (r < 1.0d0)
                atr = ptr * (1.0d0 + r * acos(r) - sin(acos(r)));
            else
                atr = ptr;
            end
        else
            %            nighttime, ptr assumed constant over nighttime
            if (supply <= 0.0d0 || ptr <= 0.0d0)
                atr = 0.0d0;
            else
                atr = min(supply, ptr);
            end
        end
        %        distribute total transpiration rate to layers
        for i = 1:nlayer
            if (flag(i) == 1)
                atrani(i) = 0.0d0;
            else
                atrani(i) = ((psiti(i) - psit) / 1000.0d0 + rt * atr) /ri(i);
                %               check for any negative transpiration losses
                if (atrani(i) < -.0000010d0)
                    negflag = 1;
                end
            end
        end
        if (nooutf == 1 && negflag == 1)
            %            find layer with most negative transpiration and omit it
            idel = 0;
            trmin = 0.0d0;
            for i = 1:nlayer
                if (atrani(i) < trmin)
                    trmin = atrani(i);
                    idel = i;
                end
            end
            flag(idel) = 1;
            %            repeat main loop with flagged layers excluded
        else
            %           done
        end
    end
end
end

%%
function  vrfli = vert(kk, kk1, ksat, ksat1, thick, thick1, psit, psit1, stone, stone1, rhowg)
% vertical flow rate
% flow rate = gradient * cond    / rhog
% mm/day   = kpa/mm   * mm/day  / kpa/mm

% Input
%     kk       hydraulic conductivity for upper layer, mm/d
%     kk1      hydraulic conductivity for lower layer, mm/d
%     ksat     saturated hydraulic conductivity of upper layer,mm/d
%     ksat1    saturated hydraulic conductivity of lower layer,mm/d
%     thick    thickness of upper layer, mm
%     thick1   thickness of lower layer, mm
%     psit     total potential of upper layer, kpa
%     psit1    total potential of lower layer, kpa
%     stone    stone volume fraction of upper layer, unitless
%     stone1   stone volume fraction of lower layer, unitless
%     rhowg    density of water times gravity acceleration, kpa/mm

% output
%     vrfli   ! vertical drainage rate from layer i, mm/d

 
    kkmean = exp((thick1 * log(kk) + thick * log(kk1)) / (thick + thick1));
%      limit kkmean to lesser saturated conductivity
    if (kkmean > ksat)
        kkmean = ksat;
    end
    if (kkmean > ksat1)
        kkmean = ksat1;
    end
    grad = (psit - psit1) / ((thick + thick1) / 2.0d0);
    vrfli = (grad * kkmean / rhowg) * (1.0d0 - (stone + stone1) / 2.0d0);
 
end

%% END MODEL





