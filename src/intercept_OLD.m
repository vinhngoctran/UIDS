%============================================================================
%   intercept.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is the Canopy interception and Snowmeil model based on 
%   Rutter–Gash canopy formulation and temperature-based model.
%   This code is developed based on Python code of Simon Dadson (2020)
%   for Artemis model.
%============================================================================

function [canopystore, snowpackstore, throughfall, snowmetl, et_soil, et_canopy, rain_int]=intercept(deltaT,rain_int,...
                                                forcing,vegheight,lai,canopystore,snowpackstore)
%% Input variables  

% rain_int                                                      % Precipitation flux [m/s]
humidity = forcing(1);                                          % Specific humidity [l]
temperature = forcing(3);                                       % Air temperature [K]    1C = 274.15K
temperature = temperature + 273.15;
windspeed = forcing(6);                                         % Wind speed [m s-1]
sradiation = forcing(4);                                        % Surface downwelling shortwave flux in air ['W m-2]                                       
lradiation = forcing(5);                                        % Surface downwelling longwave flux in air ['W m-2]       
albedo = forcing(2);                                            % Surface albedo [l]
% vegheight                                                     % Vegetation height [m]
% lai                                                           % Leaf area index (LAI) [-]
% canopystore                                                   % Canopy storage at previous time step
% snowpackstore                                                 % Storage storage at previous time step

% Constant variables:
G_sigmasb = 5.67e-8;                                              % Stefan-Boltzman constant [W m-2 K-4]
G_kappa = 0.41;                                                   % Von Karman constant [-] = [0-1]
G_lambdalh = 2.48e6;                                              % Latent heat of vaporization [J kg-1]
G_rhoa = 1.22;                                                    % Specific mass of (dry) air [kg m-3]
G_Cp = 1.01e3;                                                    % Specific heat of air [J kg-1 K-1]
G_Rv = 461.5;                                                     % Gas constant for water vapour [J kg-1 K-1]
G_pfreez = 273;                                                   % Freezing point of water [K]
G_Gam_ma = 0.066;                                                 % Psychrometric constant [kPa/K]
G_Ds = 6e-9;                                                      % Drainage rate when canopy is full [m s-1]
G_candec = 0.062;                                                 % Canopy drainage decay term [s-1]
G_F = 4.6e-8;                                                     % Degree-day melting factor [m K-1 s-1]
G_zwd = 10;                                                       % Height measurement for wind speed [m]
G_Rc = 40;                                                        % Canopy resistance (calc from LAI) [s/m] (Beven 2000 p. 76)

%%
% Compute derived parameters
% Vegetation roughness length [m] (3e-4 is for bare soil)
if vegheight == 0
    vegheight = 3e-4;
else
    vegheight = 0.1 * vegheight;
end
% zero-plane displacement [m]
d = 0.75 * vegheight;
% Temperature for snowfall [K]
t_snow = G_pfreez;
% Snowmelt threshold [K]
T_F = G_pfreez;
%   Canopy storage capacity [m] (Hough and Jones, MORECS)
C_t = 0.002 * lai;
%  throughfall coefficient [-]
if lai == 0
    phi_t = 1;
else
    phi_t = 1 - 0.5 * lai;
end
%   Compute derived meteorological variables
T_degC = temperature - 273.;
%   Scale to 2m wind speed
windspeed = windspeed * 4.87 / (log(67.8 * G_zwd - 5.42));
%   Hendrik; kPa
e_sat = 0.6108 * exp(T_degC * 17.27 / (237.3 + T_degC));
%   kg/kg to kPa; Parish and Putnam 1977)
e_abs = humidity * G_rhoa * G_Rv * temperature / 1000;
%   slope of SVP curve kPa / degC
Delta_e = 4098. * e_sat / (237.3 + T_degC)^2;
%  Aerodynamic resistance to evaporation
r_a = log((G_zwd - d) / vegheight)^2 / (G_kappa^2 * windspeed);
%   Prevent negative r_a
if r_a <= 0
    r_a = 1;
end
% Calculate upwelling shortwave "rsus" with albedo
rsus = albedo * sradiation;
% Calculate upwelling longwave "rlus" with air temperature (assuming surface temperature is air temperature)
rlus = temperature^4 * G_sigmasb;
% Total net radiation in W m-2
R_n = sradiation + lradiation - rsus - rlus;
% Partition precipitation between canopy, snow, and soil
if temperature <= t_snow
    p_snow = rain_int;
    p_soil = 0;
    p_can = 0;
else
    p_snow = 0;
    p_soil = phi_t * rain_int;
    p_can = rain_int - p_soil;
end

%% CANOPY
% Add canopy rainfall to store from previous timestep
canopy = canopystore + p_can * deltaT;
% Calculate canopy evaporation, e_can
% Vapour pressure gradient VP, (kPa)
VP = e_sat - e_abs;
if VP < 0
    VP = 0;
end
% Penman-Monteith (rc = 0; Rutter/Gash)
e_can = (Delta_e * R_n + (G_rhoa * G_Cp * VP) / r_a) / (Delta_e + G_Gam_ma) / G_lambdalh / 1000; 
% Scale canopy evaporation by canopy wetted fraction
e_can = e_can * (canopy / C_t);
% Limit e_can to available canopy water
e_can = min(canopy / deltaT, e_can);
canopy = canopy - e_can * deltaT;
% Avoid small roundoff values
if canopy < 1.e-11
    canopy = 0;
end
% Calculate canopy throughfall, q_t
if canopy < C_t
    q_t = 0;
else
    q_t = G_Ds * exp(G_candec * (canopy - C_t));
end
% Limit throughfall if it would deplete more than available
q_t = min(canopy / deltaT, q_t);
canopy = canopy - q_t * deltaT;
% Avoid small roundoff values
if canopy < 1.e-11
    canopy = 0;
end
% Ensure canopy storage doesnt exceed capacity
if canopy > C_t
    q_t = q_t + (canopy - C_t) / deltaT;
    canopy = C_t;
end

%% SNOW PACK
% Calculate snowmelt, q_m
% Degree-day melting model
q_m = G_F * max(0., temperature - T_F);
% Update snowpack add new snow to store from previous timestep
snowpack = snowpackstore + p_snow * deltaT;
% Limit snowmelt to available snowpack
q_m = min(snowpack / deltaT, q_m);
snowpack = snowpack - q_m * deltaT;
% avoid small roundoff values
if snowpack < 1.e-11
    snowpack = 0;
end

%% SURFACE
% Surface store: add new rain, snowmelt and throughfall
surface = (p_soil + q_m + q_t) * deltaT;
% Surface evaporation: Penman-Monteith
e_surf = ((Delta_e * R_n + (G_rhoa * G_Cp * VP) / r_a)/(Delta_e + G_Gam_ma * (1 + G_Rc / r_a))) / G_lambdalh / 1000;
% Subtract canopy evaporation if possible, to conserve energy
e_surf = max(e_surf - e_can, 0.);
% Limit e_surf to available water
e_surf = min(surface / deltaT, e_surf);

%% MODEL RESULTS
canopystore = canopy;                                       % Canopy store
snowpackstore = snowpack;                                   % Snowpack store
throughfall = (q_t + p_soil);                               % Throughfall
snowmetl = q_m;                                             % Snowmetl
et_soil = e_surf;                                           % Evaporation soil surface
et_canopy = e_can;                                          % Evaporation canopy
rain_int = throughfall + snowmetl - et_soil - et_canopy;    % Effective rainfall
end