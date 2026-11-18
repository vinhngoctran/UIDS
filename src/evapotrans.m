%============================================================================
%   evapotrans.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program uses an energy balance approach to calculate evapotranspiration.
%   The programe is developed based on HMRET model from:
%   Zipper, S.C. & S.P. Loheide II (2014), DOI: 10.1016/j.agrformet.2014.06.009
%============================================================================

function [act_evap,evap,diagnostics] = evapotrans(Start_Date,deltaT,t, x_cor,y_cor,n_row,x,y,s_grid,areazone,...
                forcing,lai, vegheight, tcan, albs, albv, emis, emiv,rain_int,oflow,stor,...
                G_daylightSavings, G_Zair, G_Zu, G_stefBoltz, G_Rv, G_cP, G_mW, G_mA, G_rhow, G_esA, G_esB, G_esC, G_kappa)
% First output is potential ET, retained in its historical output position.
% Second output is actual surface evaporation, limited by available water.
validateattributes(forcing,{'numeric'},{'vector','numel',6,'real','finite'});
validateattributes(deltaT,{'numeric'},{'scalar','real','finite','positive'});
validateattributes([lai vegheight rain_int oflow],{'numeric'},{'real','finite','nonnegative'});
validateattributes([albs albv emis emiv],{'numeric'},{'real','finite','>=',0,'<=',1});
validateattributes(tcan,{'numeric'},{'scalar','real','finite','>',-273.15});

%% Input variables

% Forcings
TA = forcing(1);                                            % Air temperature, degrees celsius
TA_K = TA+273.15;                                           % Air temperature, in Kelvin [degK]
SR = forcing(2);                                            % Incoming shortwave radiation, W/m2
WD = forcing(3);                                            % Wind speed, m/s
VP = forcing(4);                                            % Air vapor pressure, kPa
RH = forcing(5);                                            % Specific humidity, kg/kg (legacy column RH)
PA = forcing(6);                                            % Atmospheric pressure, kPa

% Sites
% longitude = Longitude, decimal degrees
% latitude = Latitude, decimal degrees
% albSoil = Albedo of soil, 0-1
% albVeg = Albedo of vegetation, 0-1
% emissSoil = Emissivity of soil, 0-1
% emissVeg = Emissivity of vegetation, 0-1

% Canopy
% lai = Leaf Area Index, m2/m2
% vegheight = Canopy height, m
% Tcan = Canopy surface temperature, degrees celsius
Tcan_K = tcan+273.15;                                           % Surface temperature [degK]
rhoA = (44.6*PA*273.15)/(101.3*TA_K);                           % Molar density of air, C&N eq. 3.3 [mol m-3]

% Constant variables (load common)
% global G_daylightSavings G_Zair G_Zu G_stefBoltz G_Rv G_cP G_cP G_mW G_mA G_rhow G_esA G_esB G_esC G_kappa

stime = datenum(Start_Date + seconds(deltaT) * t);
[lat,long] = coordinate(x_cor + y * s_grid - s_grid/2,y_cor + n_row * s_grid - x * s_grid + s_grid/2,areazone);
dateVector = datevec(stime);                                    % Vector of datetime
year = dateVector(1);                                           % Separate out year
julianDay = floor(stime - datenum(year,01,01))+1;               % Julian day (1 = Jan 1, etc.)
julianTime = dateVector(4)+dateVector(5)/60+dateVector(6)/3600;  % Julian time, hours (0-24)

%% Data processing, define constants
if any(~isfinite([forcing(:); lai; vegheight; tcan])) || TA_K <= 0 || Tcan_K <= 0 || PA <= 0
    error('OFM:ETInput','ET requires finite forcing and positive absolute temperature and pressure.');
end
WD = max(WD,0.1);
SR = max(SR,0);
% The legacy RH column is specific humidity [kg/kg], not relative humidity.
if VP == -9999
    if RH < 0 || RH >= 1, error('OFM:Humidity','Specific humidity must lie in [0,1).'); end
    epsilon = G_mW/G_mA;
    VP = RH*PA/(epsilon+(1-epsilon)*RH);
end
if VP < 0, error('OFM:Humidity','Vapour pressure must be nonnegative or -9999.'); end
VP = max(VP,0.01);

if vegheight < 0.01; vegheight = 0.01; end

%% Calculate Energy Balance
%% Calculate R
% Calculations for sun position from Campbell & Norman, chapter 11
sp = 279.575+0.9856*julianDay;
eqTime = (1/3600)*(-104.7*sind(sp)+596.2*sind(2*sp)+4.3*sind(3*sp)-...
    12.7*sind(4*sp)-429.3*cosd(sp)-2*cosd(2*sp)+19.3*cosd(3*sp));   % equation of time (C&N 11.4)
% coordinate() returns east-positive longitude. Use the nearest standard
% meridian when no explicit time-zone meridian is supplied by this interface.
LC = long-15*round(long/15);
solarNoon = 12+G_daylightSavings-LC/15-eqTime;
solarDec = asind(0.39785*sind(278.97+0.9856*julianDay+ ...
    1.9165*sind(356.6+0.9856*julianDay)));                          % [deg] - solar declination angle (C&N eq. 11.2)
zenith = acosd(max(-1,min(1,sind(lat)*sind(solarDec)+ ...
    cosd(lat)*cosd(solarDec)*cosd(15*(julianTime-solarNoon)))));      % [deg] - solar zenith angle (C&N eq. 11.1)

% Calculate downwelling LW radiation, following approach of Crawford &
% Duchon (1999). Includes cloudiness effects.
solConst = 1361.5;    % [W m-2] - average annual solar constant. source = wikipedia
cosZenith = max(0,cosd(zenith));
m = 35*cosZenith/sqrt(1224*cosZenith^2+1);      % calculate airmass number

% look up table for G_Tw (delta) in Tw calculation, from Smith (1966) Table 1
% Using only summer values here, but table also has spring/winter/fall
%    (potential future improvement- select values based on julianDay)
absLatitude=abs(lat);
if absLatitude < 10;
    G_Tw = 2.80;
elseif absLatitude <20;
    G_Tw = 2.70;
elseif absLatitude <30;
    G_Tw = 2.98;
elseif absLatitude <40;
    G_Tw = 2.92;
elseif absLatitude <50;
    G_Tw = 2.77;
elseif absLatitude <60;
    G_Tw = 2.67;
elseif absLatitude <70;
    G_Tw = 2.61;
elseif absLatitude <80;
    G_Tw = 2.24;
else
    G_Tw = 1.94;
end
Tdew = G_esC*log(VP/G_esA)/(G_esB - log(VP/G_esA));                 % [degC] - dewpoint temperature, from C&N eq. 3.14
atmWater = exp(0.1133 - log(G_Tw+1) + 0.0393*Tdew);                 % atmospheric precipitable water, from Crawford & Duchon

TrTpg = 1.021 - 0.084*(m*(0.000949*PA*10+0.051))^0.5;               % C&D Eq. 8 - corrections for Rayleigh scattering, absorption by permanent gases
Tw = 1 - 0.077*(atmWater*m)^0.3;                                    % C&D Eq. 9 - correction for absorption by water vapor
Ta = 0.935^m;                                                       % C&D Eq. 10 - correction for scattering by aerosols

% Calculate LWin based on cloudiness
Rso = max(0,solConst*cosZenith*TrTpg*Tw*Ta);                            % [W m-2] - clear sky shortwave irradiance
if Rso > 1
    clf = max(0,min(1,1-SR/Rso));
else
    % Shortwave observations cannot infer nighttime cloud cover. Use the
    % clear-sky emissivity term; callers with measured LW use intercsnow.
    clf = 0;
end                                 % [-] - cloudiness fraction, from 0-1. Crawford & Duchon (1999)
emissSky = max(0,min(1,(clf + (1-clf)*(1.22+0.06*sin((dateVector(2)+2)*pi/6))*(VP*10/TA_K)^(1/7))));    % [-] - emissivity based on cloud fraction, from Crawford & Duchon Eq. 20
LWin = emissSky*G_stefBoltz*(TA_K^4);                                 % [W m-2] - total incoming absorbed longwave radiation from atmosphere

% Calculate LWout based on separate vegetation & soil components
% (two-source model)
fc = 1-exp(-0.5*lai);                                               % Norman et al. (1995) Eq. 3 - fractional plant cover based on LAI
LWoutVeg = emiv*G_stefBoltz*(Tcan_K^4)*(1-exp(0.9*log(1-fc)));        % [W m-2] - outgoing LW from vegetation
LWoutSoil = emis*G_stefBoltz*(Tcan_K^4)*exp(0.9*log(1-fc));           % [W m-2] - outgoing LW from soil
LWout = LWoutVeg + LWoutSoil;                                       % [W m-2] - total outgoing longwave radiation

% Calculate total SWout as the sum of vegetation & soil components.
% (two-source model)
SWoutVeg = SR*(1-exp(0.9*log(1-fc)))*albv;                          % [W m-2] - outgoing shortwave radiation from vegetation canopy, based on amount of SW radiation reaching ground (Norman et al. (1995) Eq. 13))
SWoutSoil = SR*exp(0.9*log(1-fc))*albs;                             % [W m-2] - outgoing shortwave radiation from soil surface, based on amount of SW radiation reaching ground (Norman et al. (1995) Eq. 13))
SWout = SWoutVeg + SWoutSoil;                                       % [W m-2] - outgoing SW radiation as sum of soil and vegetation components

% Calculate net radiation budget (R)
R = SR-SWout+LWin-LWout;                                            % [W m-2] - net radiation at surface


%% Calculate Ground heat flux (G) based on amount of radiation reaching the ground
G = 0.35.*R.*exp(0.9.*log(1-fc)); % Eq. 13 - calculate G as 35% of R reaching soil (Norman et al., 1995)

%% Iterative H calculation

%Raupach (1994) z0m, d values as function of LAI, h
Cw = 2;     % empirical coefficient
Cr = 0.3;   % empirical coefficient
Cs = 0.003; % empirical coefficient
Cd1 = 7.5;  % empirical coefficient
uMax = 0.3; % empirical coefficient
subRough = log(Cw)-1+(1/Cw); % roughness-sublayer influence function

roughnessLAI=max(lai,0.01); % Bare-ground limit for roughness, not radiation cover.
u_uh = min([uMax ((Cs+Cr*(roughnessLAI/2))^0.5)]);    % ratio of u*/uh
d = vegheight*(1-(1-exp(-sqrt(Cd1*roughnessLAI)))/sqrt(Cd1*roughnessLAI)); % [m] displacement
z0m = vegheight*(1-d/vegheight)*exp(-G_kappa*(1/u_uh)-subRough);          % [m] - roughness length for momentum transfer
kB1 = 2.3;    % kB^-1 factor from Bastiaansen SEBAL paper to convert from z0m to z0h; kB1=2.3 means z0h = 0.1*z0m, which corresponds to C&N empirical equation
z0h = z0m/exp(kB1); % [m] - roughness length for heat transfer

% Solve the original sensible-heat/stability closure with a bounded root
% search. No unbounded iteration or division by a zero previous heat flux.
if G_Zu <= d+z0m || G_Zair <= d+z0h
    error('OFM:ETGeometry','Measurement heights must exceed displacement and roughness heights.');
end
[H,stability,iterations] = sensible_heat(tcan-TA,WD,rhoA,TA_K,d,z0m,z0h, ...
    G_Zu,G_Zair,G_kappa,G_cP,G_mA);

%% Calculate ET rate!
ET = R-G-H;                                                         % [W m-2] - ET rate as residual of energy budget
lamda = (2.495-(2.36e-3)*tcan)*G_mW*1000000;                        % [J mol-1] - latent heat of vaporization of water, dependent on temp [degC]. Formula B-21 from Dingman 'Physical Hydrology'
act_evap = (G_mW/G_rhow)*(60*60)*1000*ET/lamda;                     % [mm hr-1] - evaporation rate
act_evap = max(act_evap, 0);                                        % set equal to 0 if it calculates something below 0 (seems to happen occasionally in areas with very low veg cover)
act_evap = act_evap / 1000/ 3600;                                   % Convert mm hr-1 to m s-1
evap = min(act_evap,max(0,rain_int+max(0,oflow)/deltaT));
diagnostics = struct('NetRadiation',R,'GroundHeat',G,'SensibleHeat',H, ...
    'Latitude',lat,'Longitude',long,'VapourPressure',VP,'Stability',stability, ...
    'Iterations',iterations,'PotentialET',act_evap,'SurfaceET',evap);         % Evapotranspiration [m s-1]
end

function [H,zeta,iterations] = sensible_heat(dT,wind,rho,T,d,z0m,z0h,Zu,Za,kappa,cP,mA)
if abs(dT) < 1e-12, H=0; zeta=0; iterations=0; return; end
% The root has the opposite sign from surface-minus-air temperature.
lo=0; hi=sign(-dT);
f0=residual(0); f1=residual(hi);
for k=1:60
    if sign(f0) ~= sign(f1), break; end
    hi=2*hi; f1=residual(hi);
end
if sign(f0) == sign(f1)
    error('OFM:ETConvergence','Could not bracket atmospheric stability. Check forcing and canopy geometry.');
end
if hi < lo, temp=lo; lo=hi; hi=temp; end
flo=residual(lo);
for iterations=1:100
    zeta=(lo+hi)/2;
    f=residual(zeta);
    if abs(f) <= 1e-10*max(1,abs(zeta)), break; end
    if sign(f)==sign(flo), lo=zeta; flo=f; else, hi=zeta; end
end
[~,H]=residual(zeta);
if ~isfinite(H) || ~isreal(H)
    error('OFM:ETConvergence','Invalid sensible heat flux.');
end
    function [f,heat] = residual(z)
        if z > 0
            dm=6*log1p(z); dh=dm;
        else
            dh=-2*log((1+sqrt(1-16*z))/2); dm=0.6*dh;
        end
        lm=max(1e-6,log((Zu-d)/z0m)+dm);
        lh=max(1e-6,log((Za-d)/z0h)+dh);
        ustar=wind*kappa/lm;
        rAir=1.4*lm*lh/(kappa^2*wind*rho);
        rExcess=mA*log(z0m/z0h)/(rho*kappa*ustar);
        heat=cP*dT/(rAir+rExcess);
        f=z+kappa*9.8*Zu*heat/(rho*cP*T*ustar^3);
    end
end
