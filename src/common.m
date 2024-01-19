%============================================================================
%   common.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program concludes common variables for modules of PCUW
%============================================================================



%% Evapotranspiration Model
global G_daylightSavings G_Zair G_Zu G_stefBoltz G_Rv G_cP G_mW G_mA G_rhow G_esA G_esB G_esC G_kappa

%% Canopy interception and snowmelt models
global G_lambdalh G_rhoa G_Cp G_pfreez G_Gam_ma G_Ds G_candec G_F G_zwd G_Rc

%% Common variables
G_daylightSavings = 1;                                            % 1 if DST (summer), 0 if not
G_Zair = 3.66;                                                    % Height of air temperature measurements [m]
G_Zu = 9.144;                                                     % Height of wind speed measurements [m]
G_stefBoltz = 5.67e-8;                                            % Stefan-Boltzmann Constant [W m-2 K-4]
G_Rv = 461.5;                                                     % Gas constant for water vapour [J kg-1 K-1]
G_cP = 29.3;                                                      % Specific heat of air, from Campbell & Normal Table A1 [J mol-1 degC-1]
G_mW = 0.018;                                                     % Molecular mass of water, from C&N table A1 [kg mol-1]
G_mA = 0.029;                                                     % Molecular mass of air, from C&N table A1 [kg mol-1]
G_rhow = 1000;                                                    % Density of water, from C&N table A2 [kg m-3]
G_esA = 0.611;                                                    % Constant 'a' for es equations, C&N pg 41 [kPa]
G_esB = 17.502;                                                   % Constant 'b' for es equations, C&N pg 41 [-]
G_esC = 240.97;                                                   % Constant 'c' for es equations, C&N pg 41 [degC] -
G_kappa = 0.41;                                                   % von Karman's constant [-]
G_lambdalh = 2.48e6;                                              % Latent heat of vaporization [J kg-1]
G_rhoa = 1.22;                                                    % Specific mass of (dry) air [kg m-3]
G_Cp = 1.01e3;                                                    % Specific heat of air [J kg-1 K-1]
G_pfreez = 273;                                                   % Freezing point of water [K]
G_Gam_ma = 0.066;                                                 % Psychrometric constant [kPa/K]
G_Ds = 6e-9;                                                      % Drainage rate when canopy is full [m s-1]
G_candec = 0.062;                                                 % Canopy drainage decay term [s-1]
G_F = 4.6e-8;                                                     % Degree-day melting factor [m K-1 s-1]
G_zwd = 10;                                                       % Height measurement for wind speed [m]
G_Rc = 40;                                                        % Canopy resistance (calc from LAI) [s/m] (Beven 2000 p. 76)



