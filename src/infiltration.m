%============================================================================
%   infiltration.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate infiltration rate by Green-Amplt method follows 
%   the method of Chu (1978), Mein and Larson (1971, 1973) 
%   and Skaggs and Kaheel (1982) and based on the fortran code of 
%   J.E. Parsons and R. Mu?oz-Carpena 
%   (provided in: http://abe.ufl.edu/carpena/software/wingampt.shtml)
%============================================================================

function [bf,inf_rate,stor,ro,ipond,fp,tp,tpp,v_infil,of_rate]=...
    infiltration(bf,inf_rate,stor,ro,ipond,fp,tp,tpp,x,y,deltaT,soil,v_infil,pa_infil,oflow,of_rate,imp)

%% Variables
% bf                                % Cumulative precipitation for infiltration (cm)
% inf_rate                          % Infiltration rate (cm/h)
% stor                              % Soil moisture [cm]
% ro                                % Runoff (cm)
% ipond                             % Ponded conditions
% fp                                % Infiltration rate at ponding (cm/h)
% tp                                % Time to ponding (h)
% tpp                               % tpp = time required to infiltrate Fp if the system had started in ponded conditions (h)

%%
% Soil type
iinf = soil(x,y);
%soil variables: infitration parameters
Ks = pa_infil(iinf,1);                                                              % vertical saturated conductivity (cm/h)
Sw = pa_infil(iinf,2);                                                              % average suction across the wetting front (cm))
Wh = pa_infil(iinf,3);                                                              % water content, theta at saturation (cm^3/cm^3)
Whin = pa_infil(iinf,4);                                                            % initial water content, theta (cm^3/cm^3) at the start
Stm = pa_infil(iinf,5);                                                             % maximum surface storage s(cm)

% Previous infiltration rate
inf_rate = inf_rate * 100 * 3600;
% Compute rainfall
rain_in = (of_rate * (1 - imp(x,y))) * 100 * 3600  + oflow * 100;          % Convert m/s to cm/hour
bm = Wh - Whin;
deltaT = deltaT/3600;

%   Not ponded at the start of the period, then
if ipond <= 0
    
    %  bfp: potential infiltration rate
    if fp > 0 && rain_in >= fp
        tp = 0;
        bfp = Sw * bm*Ks/(fp-Ks);
    elseif rain_in > Ks
        bfp = Sw * bm /((rain_in/Ks)-1.0);
        tp = bfp / rain_in;
    else
        bfp = 0;
        tp = 9999;
    end
    
    % tpp: time shift   
    tpp = (bfp - Sw * bm * log(1.0 + bfp / (Sw * bm))) / Ks;   
    
    if tp > deltaT      
        % no ponding during this period
        tp = 9999;
        tpp = 9999;        
        bf = bf + rain_in * deltaT;
        fp = Ks + (Ks*bm*Sw/bf);
        inf_rate = rain_in;
        stor = 0.0;
        ipond=  0;
    else     
        % ponding at tp
        bf = bf + rain_in * deltaT;
        fp = Ks + (Ks*bm*Sw/bf);
        inf_rate = rain_in;
        stor = 0.0;
        
        %* ponded at deltaT
        [bf, inf_rate, stor, ro] =  pndinf(deltaT,rain_in,tp,tpp,Ks, Sw, bm, Stm,bf, stor, ro);
        ipond=1;   
    end
else    % if the period starts with ponding, then
    
    %  find time to infiltrate fnp
    if rain_in < inf_rate
        amtinf = bf + stor;        
        tnp = newtnp(deltaT,tp,tpp,rain_in,amtinf,Ks, Sw, bm);    
    else
        
        % will not loose ponding, set tnp>deltaT
        tnp = deltaT+1.0;
    end
    
    if tnp > deltaT       
        % ponding for whole period
        [bf, inf_rate, stor, ro] = pndinf(deltaT,rain_in,tp,tpp,Ks, Sw, bm, Stm,bf, stor, ro);
        ipond = 1;
    else
        % ponding ends at tnp: ponded portion
        [bf, inf_rate, stor, ro] = pndinf(deltaT,rain_in,tp,tpp,Ks, Sw, bm, Stm,bf, stor, ro);
        
        % no pond portion
        bf = bf + rain_in * deltaT;
        fp = Ks + (Ks * bm *Sw / bf);
        inf_rate = rain_in;
        stor = 0.0;    
        ipond = 0;
    end
end
% Convert infiltration rate from cm/h to m/s
inf_rate = inf_rate / 100 / 3600; 
% Update infiltration depth
v_infil = v_infil + inf_rate * deltaT;
of_rate = of_rate - inf_rate ;
end

%%
function tnp = newtnp(deltaT,tp,tpp,rain_in,amtinf,Ks, Sw, bm)

accpt = 0.1e-10;
iter = 0;
tnp = deltaT;
bftry = (tnp) * rain_in + amtinf;
arglog = 1.0 + bftry/(bm*Sw);
hnp = bftry - bm*Sw*log(arglog) -  Ks*(tnp -tp +tpp);
error = 1;

while error > accpt
    iter = iter+1;
    tnpold = tnp;
    if iter > 1000
        break;
    end    
    dhnp= rain_in - rain_in*(1.0/arglog) - Ks;  
    tnp = tnpold - hnp/dhnp;   
    % if tnp is negative  - fix added jep
    if tnp < 0
        tnp = 0;
    end
    bftry = (tnp)*rain_in + amtinf;
    arglog = 1.0 + bftry/(bm*Sw);    
    hnp = bftry - bm*Sw*log(arglog) - Ks*(tnp -tp +tpp);
    error = abs(hnp);
end
end

%%
function [bf, inf_rate, stor, ro] =  pndinf(deltaT,rain_in,tp,tpp,Ks, Sw, bm, Stm,bf, stor, ro)

waterdepth = rain_in * deltaT + stor;
% make a guess for bigf
bbf = bf+waterdepth;
bbf = sschu(deltaT,tp,tpp,bbf,Ks, Sw, bm);
delinf = bbf - bf;
bf = bbf;
inf_rate = Ks + (Ks*bm*Sw/bf);
if waterdepth > delinf
    stor = waterdepth - delinf;
    if stor > Stm
        ro = ro + (stor-Stm);
        stor = Stm;
    end
else
    stor = 0.0;
end
end

%%
function bff = sschu(deltaT, tp, tpp, bff, Ks, Sw, bm)
% use Newton's method on chu's equation
iter = 0;
accpt = 0.1e-04; % threshold
hh = bff - bm * Sw * log (1.0 + (bff/(bm*Sw))) - Ks * (deltaT-tp+tpp);
error = 1;
while error > accpt
    iter = iter+1;
    if iter > 200
        break;
    end
    dhdf = 1.0 - ((bm * Sw) /(bm*Sw + bff));
    bff = bff - hh/dhdf ;
    hh = bff - bm * Sw * log (1.0 + (bff/(bm*Sw))) - Ks * (deltaT-tp+tpp);
    error = abs(hh);
end
end
