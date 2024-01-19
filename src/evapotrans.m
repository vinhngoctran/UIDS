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

function [act_evap evap] = evapotrans(Start_Date,deltaT,t, x_cor,y_cor,n_row,x,y,s_grid,areazone,...
                forcing,lai, vegheight, tcan, albs, albv, emis, emiv,rain_int,oflow,stor,...
                G_daylightSavings, G_Zair, G_Zu, G_stefBoltz, G_Rv, G_cP, G_mW, G_mA, G_rhow, G_esA, G_esB, G_esC, G_kappa)

%% Input variables

% Forcings
TA = forcing(1);                                            % Air temperature, degrees celsius
TA_K = TA+273.16;                                           % Air temperature, in Kelvin [degK]
SR = forcing(2);                                            % Incoming shortwave radiation, W/m2
WD = forcing(3);                                            % Wind speed, m/s
VP = forcing(4);                                            % Air vapor pressure, kPa
RH = forcing(5);                                            % Relative humidity, kg/kg
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
Tcan_K = tcan+273.16;                                           % Surface temperature [degK]
rhoA = (44.6*PA*273.15)/(101.3*TA_K);                           % Molar density of air, C&N eq. 3.3 [mol m-3]

% Constant variables (load common)
% global G_daylightSavings G_Zair G_Zu G_stefBoltz G_Rv G_cP G_cP G_mW G_mA G_rhow G_esA G_esB G_esC G_kappa

stime = datenum(Start_Date + seconds(deltaT) * t);
[long lat] = coordinate(x_cor + x * s_grid - s_grid/2,y_cor + n_row * s_grid - y * s_grid + s_grid/2,areazone);
dateVector = datevec(stime);                                    % Vector of datetime
year = dateVector(1);                                           % Separate out year
julianDay = floor(stime - datenum(year,01,01))+1;               % Julian day (1 = Jan 1, etc.)
julianTime = dateVector(4)+(dateVector(5)/60);                  % Julian time, in hours (0-24)

%% Data processing, define constants
% Check if climate variables <= 0
if WD < 0.1
    WD = 0.1; 
end
if VP < 0.01
    VP = 0.01;
elseif VP == -9999                                                  % calculated from Tair and relative humidity
    % Hendrik; kPa
    e_sat = 0.6108 * exp(TA * 17.27 / (237.3 + TA));
    % kg/kg to kPa; Parish and Putnam 1977)
    e_abs = RH * rhoA * G_mA * G_Rv * TA / 1000;
    % Vapour pressure gradient VP, (kPa)
    VP = e_sat - e_abs;
end

if lai < 0.01; lai = 0.01; end
if vegheight < 0.01; vegheight = 0.01; end

%% Calculate Energy Balance
%% Calculate R
% Calculations for sun position from Campbell & Norman, chapter 11
sp = 279.575+0.9856*julianDay;
eqTime = (1/3600)*(-104.7*sind(sp)+596.2*sind(2*sp)+4.3*sind(3*sp)-...
    12.7*sind(4*sp)-429.3*cosd(sp)-2*cosd(2*sp)+19.3*cosd(3*sp));   % equation of time (C&N 11.4)
while long > 7.5;
    long = long - 15;                                               % figure out # of degrees west of standard meridian (if negative, that means you're east!)
end
LC = -long/15;                                                      % calculate longitude correction
solarNoon = 12 + G_daylightSavings - LC - eqTime;                     % [hrs] - calculate solar noon time
solarDec = asind(0.39785*sind(278.97+0.9856*julianDay+ ...
    1.9165*sind(356.6+0.9856*julianDay)));                          % [deg] - solar declination angle (C&N eq. 11.2)
zenith = acosd(sind(lat)*sind(solarDec)+ ...
    cosd(lat)*cosd(solarDec)*cosd(15*(julianTime-solarNoon)));      % [deg] - solar zenith angle (C&N eq. 11.1)

% Calculate downwelling LW radiation, following approach of Crawford &
% Duchon (1999). Includes cloudiness effects.
solConst = 1361.5;    % [W m-2] - average annual solar constant. source = wikipedia
m = 35*cosd(zenith)*(1224*cosd(zenith)*cosd(zenith)+1)^(-0.5);      % calculate airmass number

% look up table for G_Tw (delta) in Tw calculation, from Smith (1966) Table 1
% Using only summer values here, but table also has spring/winter/fall
%    (potential future improvement- select values based on julianDay)
if lat < 10;
    G_Tw = 2.80;
elseif lat <20;
    G_Tw = 2.70;
elseif lat <30;
    G_Tw = 2.98;
elseif lat <40;
    G_Tw = 2.92;
elseif lat <50;
    G_Tw = 2.77;
elseif lat <60;
    G_Tw = 2.67;
elseif lat <70;
    G_Tw = 2.61;
elseif lat <80;
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
Rso = solConst*cosd(zenith)*TrTpg*Tw*Ta;                            % [W m-2] - clear sky shortwave irradiance
clf = max([0 min([(1-SR/Rso) 1])]);                                 % [-] - cloudiness fraction, from 0-1. Crawford & Duchon (1999)
emissSky = (clf + (1-clf)*(1.22+0.06*sin((dateVector(2)+2)*pi/6))*(VP*10/TA_K)^(1/7));    % [-] - emissivity based on cloud fraction, from Crawford & Duchon Eq. 20
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

u_uh = min([uMax ((Cs+Cr*(lai/2))^0.5)]);    % ratio of u*/uh
d = vegheight*(1-(1-exp(-sqrt(Cd1*lai)))/sqrt(Cd1*lai));            % [m] - zero-plane displacement height
z0m = vegheight*(1-d/vegheight)*exp(-G_kappa*(1/u_uh)-subRough);          % [m] - roughness length for momentum transfer
kB1 = 2.3;    % kB^-1 factor from Bastiaansen SEBAL paper to convert from z0m to z0h; kB1=2.3 means z0h = 0.1*z0m, which corresponds to C&N empirical equation
z0h = z0m/exp(kB1); % [m] - roughness length for heat transfer

%%%%%%%%%%%%%% Begin Iteration - positive stability %%%%%%%%%%%%%%%%%%%
% Iterative solution to H, rH, etc. starting from very low values of H and positive stability
zeta = 0.5;                                                         % initial stability factor for diabatic correction (zeta from C&N sec 7.5) - >0 when surface cooler than air
Hiter = 0.5;                                                        % Hiter is placeholder for H during iterative process
changePerc = 0.5;                                                   % arbitrary starting value;
i = 0;                                                              % starting i for iterations
while abs(changePerc)>0.001                                         % set convergence criteria here, in percent
    i = i+1;                                                        % advance iteration number by one
    % calculate diabatic correction factors based on zeta. from C&N
    if zeta > 0                                                     % stable flow
        %   C&N equation 7.27
        diabM = 6*log(1+zeta);                                      % diabatic correction for momentum transfer
        diabH = diabM;                                              % diabatic correction for heat transfer
    else
        %    C&N equation 7.26
        diabH = -2*log((1+(1-16*zeta).^0.5)/2);
        diabM = 0.6*diabH;
    end
    
    %  calculate u*, gHa based on diabatic correction factors. from C&N
    uStar = WD*G_kappa/(log((G_Zu-d)/z0m)+diabM); % [m] - friction velocity C&N eq. 7.24
    rHa = 1/((G_kappa^2)*WD*rhoA/(1.4*((log((G_Zu-d)/z0m)+diabM)*(log((G_Zair-d)/z0h)+diabH)))); % [m2 s mol-1]
    
    % by including rExcess, we can ignore the difference between Tsurf and Taero
    rExcess = G_mA*log(z0m/z0h)/(rhoA*G_kappa*uStar);                       % [m2 s mol-1] - excess resistance, from Norman & Becker (1995)
    
    rHtot = rHa+rExcess; % [m2 s mol-1] - total resistance
    
    % Calculate H and zeta for next iteration
    Hiter(i,1) = G_cP*(tcan-TA)/rHtot;                                % sensible heat flux
    zeta = -G_kappa*9.8*G_Zu*Hiter(i)/(rhoA*G_cP*TA_K*(uStar.^3));            % updated zeta value
    
    % calculate the percent change between iterations
    if i>2
        changePerc = (Hiter(i)-Hiter(i-1))/Hiter(i-1);
    else
        changePerc = 0.5;                                           % Arbitrary initial value for changePerc on first iteration
    end
    
    if i == 10000                                                   % non-convergence scenario. occasionally happens when u is very very very low (or negative).
%         error('10000 iterations, will not converge');
    end
end

H_low = Hiter(length(Hiter));                                       % [W m-2] - converged H for low starting values

if H_low == 0;
    H_low = 0.02;                                                   % set to 0.02 so you don't divide by 0 later on
end

clear Hiter changePerc zeta

%%%%%%%%%%%%%% Begin Iteration - negative stability %%%%%%%%%%%%%%%%%%%
% Repeat iterative solution to H & rH, starting from very high values of H and negative stability
zeta(1) = -0.5;  % initial stability factor for diabatic correction (zeta from C&N sec 7.5) - >0 when surface cooler than air
Hiter(1) = 500; % Hiter is placeholder for H during iterative process
changePerc = 0.5; % arbitrary starting value;
i = 0;  % starting i for iterations
while abs(changePerc)>0.001                                         % set convergence criteria here, in percent
    i = i+1;                                                        % advance iteration number by one
    % calculate diabatic correction factors based on zeta. from C&N
    if zeta > 0                                                     % stable flow
        %   C&N equation 7.27
        diabM = 6*log(1+zeta);                                      % diabatic correction for momentum transfer
        diabH = diabM;                                              % diabatic correction for heat transfer
    else
        %    C&N equation 7.26
        diabH = -2*log((1+(1-16*zeta).^0.5)/2);
        diabM = 0.6*diabH;
    end
    
    % calculate u*, gHa based on diabatic correction factors. from C&N
    uStar = WD*G_kappa/(log((G_Zu-d)/z0m)+diabM);                           % [m] - friction velocity C&N eq. 7.24
    rHa = 1/((G_kappa^2)*WD*rhoA/(1.4*((log((G_Zu-d)/z0m)+diabM)*(log((G_Zair-d)/z0h)+diabH)))); % [m2 s mol-1]
    
    % by including rExcess, we can ignore the difference between Tsurf and Taero
    rExcess = G_mA*log(z0m/z0h)/(rhoA*G_kappa*uStar);                       % [m2 s mol-1] - excess resistance, from Norman & Becker (1995)
    
    rHtot = rHa+rExcess;                                            % [m2 s mol-1] - total resistance
    
    % calculate H and zeta for next iteration
    Hiter(i,1) = G_cP*(tcan-TA)/rHtot;                                % sensible heat flux
    zeta = -G_kappa*9.8*G_Zu*Hiter(i)/(rhoA*G_cP*TA_K*(uStar.^3));            % updated zeta value
    if i>2
        changePerc = (Hiter(i)-Hiter(i-1))/Hiter(i-1);
    else
        changePerc = 0.5;
    end
    
    if i == 10000
%         error('10000 iterations, will not converge');
    end
end

H_high = Hiter(length(Hiter));                                      % [W m-2] - converged H for high starting values

if H_high == 0
    H_high = 0.02;                                                  % set to 0.02 so you don't divide by 0
end

clear Hiter changePerc zeta

if abs((H_low-H_high))/H_low <= 0.01                                % if H_low and H_high are within 1% of each other, you've converged on a universal solution!
    H = (H_low+H_high)/2;                                           % [W m-2] - take the mean of your two H calculations for your sensible heat flux
else                                                                % if both H_low and H_high converge, but the numbers aren't close to each other, something weird happening - probably imaginary numbers in your solution for some reason (negative value somewhere)
    error('H_low and H_high are too far apart!')
end

%% Calculate ET rate!
ET = R-G-H;                                                         % [W m-2] - ET rate as residual of energy budget
lamda = (2.495-(2.36e-3)*tcan)*G_mW*1000000;                        % [J mol-1] - latent heat of vaporization of water, dependent on temp [degC]. Formula B-21 from Dingman 'Physical Hydrology'
act_evap = (G_mW/G_rhow)*(60*60)*1000*ET/lamda;                     % [mm hr-1] - evaporation rate
act_evap = max(act_evap, 0);                                        % set equal to 0 if it calculates something below 0 (seems to happen occasionally in areas with very low veg cover)
act_evap = act_evap / 1000/ 3600;                                   % Convert mm hr-1 to m s-1
evap = min(act_evap,rain_int + (stor/100 + oflow) /deltaT);         % Evapotranspiration [m s-1]
end