%============================================================================
%   ofm_urban.m
%
%   Project:    Overland Flow-Urban Flood Model (OFM-Urban)
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Main program of OFM-Urban.
%============================================================================

function Results = ofm_urban(Filename, options)
% Results contains the final snapshot. Each step is saved independently.
if nargin < 2, options = struct(); end
if ~isfield(options,'Plot'), options.Plot = false; end
clc
display(['================== Overland Flow-Urban Flood Model (OFM-Urban) ====================']);
display(' ');
display(['                         © Vinh Ngoc Tran, 2021 ']);
display(['                       HYDROLAB, University of Ulsan ']);
display(' ');
display(['===================================================================================']);
display(' ');


%% #################################################################################################################################################
%
%                                                          INPUT GENERALIZATION
%
% ##################################################################################################################################################

% Read input filenames
display(['Data reading . . .']);
current_folder = pwd;
Example_file = Filename;

opts = delimitedTextImportOptions("NumVariables", 9);
opts.DataLines = [1, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["Project", "TitleNotes", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
opts = setvaropts(opts, ["Project", "TitleNotes", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Project", "TitleNotes", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"], "EmptyFieldRule", "auto");
InP_FILE = readmatrix(Example_file, opts);
clear opts
InP_FILE = fillmissing(InP_FILE,'constant',"");
% Configuration files may use Windows separators on any platform.
folderRow = find(InP_FILE(:,1) == "FOLDER");
InP_FILE(folderRow,2) = string(strrep(char(InP_FILE(folderRow,2)),char(92),filesep));
restore_folder = onCleanup(@() cd(current_folder));
clear opts Example_file


% Read data
% Read Time information ---------------------------------------------------
Start_Date = datetime(str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),2)),...     % Start date of simulation
    str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),3)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),4)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),5)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),6)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["START_TIME"]),7)));
% Start date of simulation
End_Date = datetime(str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),2)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),3)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),4)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),5)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),6)),...
    str2num(InP_FILE(find(InP_FILE(:,1) == ["END_TIME"]),7)));

deltaT = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTA_T"]),2));                        % The computational time step (in second)
DELTA_T1D = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTA_T1D"]),2));                   % The computational time step for 1D sewer model (in second)
DeltaT_Time = string(duration(0,0,deltaT),'hh:mm:ss');
t_dur = seconds(End_Date-Start_Date); % Duration [seconds]
validateattributes(deltaT,{'numeric'},{'scalar','positive','finite'});
validateattributes(t_dur,{'numeric'},{'scalar','positive','finite'});
NT_RunOF = round(t_dur/deltaT);
if abs(NT_RunOF*deltaT-t_dur) > 1e-8*max(1,t_dur)
    error('OFM:TimeStep','Simulation duration must be a whole number of DELTA_T steps.');
end
NT_WAR = str2num(InP_FILE(find(InP_FILE(:,1) == ["NT_WAR"]),2));                         % Reporting interval in seconds (step snapshots are always written).
validateattributes(NT_WAR,{'numeric'},{'scalar','positive','finite'});
NT_WAR = max(1,round(NT_WAR/deltaT));
TSP2D = str2num(InP_FILE(find(InP_FILE(:,1) == ["TSTP2D"]),2));
if ~ismember(TSP2D,[1 2]), error('OFM:Scheme','TSTP2D must be 1 or 2.'); end
% TSTP2D: 1 legacy sequential update, 2 adaptive explicit routing.

% Read model option -------------------------------------------------------
modchn = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODCHN"]),2));                         % Checking for channel
moddrain = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODDRA"]),2));                       % Checking for manhole (sewer system)
infsc = str2num(InP_FILE(find(InP_FILE(:,1) == ["INFSC"]),2));                           % Checking for free outlet (sewer system)
outsc = str2num(InP_FILE(find(InP_FILE(:,1) == ["OUTSC"]),2));                           % Checking for free outlet (sewer system or pump)
out_free = str2num(InP_FILE(find(InP_FILE(:,1) == ["OUT_FREE"]),2));                     % Checking for free outlet (sewer system)

x_coord = zeros(0,1); y_coord = zeros(0,1);
n_chout = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_CHOUT"]),2));                       % Cell number in row for outlet
for i=1:n_chout
    x_coord(i,1) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,2));              % Cell number in row for outlet
    y_coord(i,1) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,3));              % Cell number in colum for outlet
    slopchn(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,4));              % Bed slope of channel at outlet
    Qmax_out(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,5));             % Maximum flow at outlet (optional or 0)
    W_out(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,6));                % Width of channel at outlet
    D_out(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,7));                % Depth of channel at outlet
    Man_out(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,8));              % Manning's coefficients of channel at outlet
    slop_out(i) = str2num(InP_FILE(find(InP_FILE(:,1) == ["CHN_INF"])+i,9));             % Slope of cell at outlet   
end

modint = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODINT"]),2));                         % Canopy interception and snowmelt model (1: used; 0: none)
modets = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODETS"]),2));                         % Energy balance-based Evapotranspiration model (1: used; 0: none)
modswf = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODSWF"]),2));                         % Soil water fluxes with vegetation (1: used; 0: none)

modinf = str2num(InP_FILE(find(InP_FILE(:,1) == ["MODINF"]),2));                         % Infiltration model (1: used; 0: none)
n_soil = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_SOIL"]),2));                         % The number of soil stypes
if modinf == 1
    for i = 1:n_soil
        pa_infil(i,:) =  [str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_INFIL"])+i,2)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_INFIL"])+i,3)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_INFIL"])+i,4)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_INFIL"])+i,5)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_INFIL"])+i,6))];                 % Infiltration parameters (5)
        
    end
    
    n_runoff = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_SUBB"]),2));                   % The number of sub-catchment with different parameters
    for i = 1:n_runoff
        pa_runoff(i,:) =  [str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,2)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,3)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,4)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,5)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,6)),...
            str2num(InP_FILE(find(InP_FILE(:,1) == ["PA_ROFF"])+i,7))];                  % Runoff parameters (6)        
    end   
end
clearvars i;

n_boun = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_BOUN"]),2));                         % The number of soft boundary conditions
if n_boun ~= 0
    deltaT_bound = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTAT_BOUN"]),2));          % Time interval for boundary condition (in second)
    boundary_dat = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},...             % Read boundary data
                                InP_FILE{find(InP_FILE(:,1)==["BOUN"]),2}]);
    % Boundary data stay at their native interval; sampled during the run.
end

if modswf == 1 && (modinf ~= 1 || (modint ~= 1 && modets ~= 1))
    error('OFM:SoilWaterOptions','MODSWF requires MODINF and either MODINT or MODETS.');
end
if modswf ==1                                                                            % Load inputs for Soil water fluxes model
    modsfl_inf = modswfip([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["SFLX"]),2}]);
end
n_chn = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_CHN"]),2));                           % The number of channel links
max_chn = str2num(InP_FILE(find(InP_FILE(:,1) == ["MAX_CHN"]),2));                       % Maximum number of channel elements

id_rain = str2num(InP_FILE(find(InP_FILE(:,1) == ["ID_RAIN"]),2));                       % Rainfall index ('1' for raingauge OR '0' for gridded rainfall)
n_gauge = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_GAUGE"]),2));                       % The number of raingauges
deltaT_rain = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTAT_RAIN"]),2));               % Time interval for rainfall (in second)
if id_rain == 0
    rgauge_name = 'NaN';
    x_rain = 0;
    y_rain = 0;
else
    for ig = 1:n_gauge
        rgauge_name{ig} = InP_FILE{find(InP_FILE(:,1) == ["GAUGE_ID"])+ig,1};            % Rain gauge name
        x_rain(ig,1) = str2num(InP_FILE(find(InP_FILE(:,1) == ["GAUGE_ID"])+ig,2));      % Longtitude of raingauges
        y_rain(ig,1) = str2num(InP_FILE(find(InP_FILE(:,1) == ["GAUGE_ID"])+ig,3));      % Latitude of raingauges
    end
end
clearvars ig;

id_flow = str2num(InP_FILE(find(InP_FILE(:,1) == ["ID_FLOW"]),2));                       % Flow station check (0 or 1)
n_flow = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_FLOW"]),2));                         % The number of flow station
if id_flow == 1
    for nf = 1:n_flow
        locat_flow(nf,:) = [str2num(InP_FILE(find(InP_FILE(:,1) == ...                   % Information of flow station
            ["FLOW_ID"])+nf,2)),str2num(InP_FILE(find(InP_FILE(:,1) == ["FLOW_ID"])+nf,3))];
    end
else
    locat_flow = [];
end
clearvars nf;

n_node = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_NODES"]),2));                       % Flow station check (0 or 1)
n_link = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_LINKS"]),2));                       % Flow station check (0 or 1)
n_out = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_OUTS"]),2));                         % Flow station check (0 or 1)

% Read channel width and link information
if modchn == 1
    CHN_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["CHNW"]),2}]);
    if n_chn > 0
        pa_chn = CHN_DAT(1:n_chn,1:4);                                                   % Channel parameters: #, width, depth, and Manning's coefficients
    end
    
    LINK_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["LINK"]),2}]);
    for i=1:n_chn
        for j=1:2
            ad_chn(i,1:max_chn+2,j) = ...
                [LINK_DAT(2+(i-1)*2+j,1:max_chn), -1, -1];                                         % Channel element addreses
        end
    end
    clearvars CHN_DAT LINK_DAT i j;
end

% Read heading
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["ncols", "VarName2","VarName3"];
opts.VariableTypes = ["string", "string", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
opts = setvaropts(opts, ["ncols", "VarName2","VarName3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ncols", "VarName2","VarName3"], "EmptyFieldRule", "auto");
HEAD_DAT = readmatrix...
    ([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},...
    InP_FILE{find(InP_FILE(:,1)==["HEAD"]),2}], opts);                                   % Heading of shape area, elevation, land use, and soil type
areazone = [HEAD_DAT{1,2},' ',HEAD_DAT{1,3}];                                            % UTM coordinate zones
n_col = str2num(HEAD_DAT(2,2));                                                          % The number of colums
n_row = str2num(HEAD_DAT(3,2));                                                          % The number of rows
x_cor = str2num(HEAD_DAT(4,2));                                                          % Corner X of area
y_cor = str2num(HEAD_DAT(5,2));                                                          % Corner Y of area
s_grid = str2num(HEAD_DAT(6,2));                                                         % Grid resolution
nodata = str2num(HEAD_DAT(7,2));                                                         % No data
head_plot = [0 -s_grid;s_grid 0;x_cor y_cor];                                            % Heading of data for result visualization
clearvars HEAD_DAT opts;

% Read shape of area
% 0:        Land
% 1:        Channel
% 2:        Outlet
% 3:        Critial boundary
% 4:        Soft boundary
% 5:        Free boundary (zero boundary)
% 999:      Building
% -9999:    No data
SHAPE_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["SHAPE"]),2}]);
shp = SHAPE_DAT(1:n_row,1:n_col);
clearvars SHAPE_DAT;

if modchn == 1 && n_chn > 0        % Define channel
    for i = 1:n_chn
        for j =1:max_chn
            x = ad_chn(i,j,1);
            y = ad_chn(i,j,2);
            if x > 0 && y > 0
                shp(x,y) = 1;
            end
        end
    end
end
clearvars i j x y;

% Read elevation (DEM)
DEM_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["DEMM"]),2}]);
dem = DEM_DAT(1:n_row,1:n_col);                                                           % Elevation for watershed
clearvars DEM_DAT;
% Compute watershed area
W_are = sum(sum(shp>-9999))*s_grid^2/1000000;                                             % Watershed area [km2]
n_cell = sum(sum(shp>-9999));                                                             % The number of grid cell in watershed area

% Read boundary shape
if n_boun ~= 0
    BOUND_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},...                 % Read boundary data
                                InP_FILE{find(InP_FILE(:,1)==["BSHA"]),2}]);
    bshape = BOUND_DAT(1:n_row,1:n_col); 
    clearvars BOUND_DAT;
end
prescribed = shp == 4 | shp == 6;
if any(prescribed(:))
    if n_boun == 0
        error('OFM:BoundaryInput','Stage and flow boundary cells require N_BOUN > 0.');
    end
    ids = bshape(prescribed);
    if any(ids < 1 | ids ~= round(ids) | ids > size(boundary_dat,2))
        error('OFM:BoundaryInput','Boundary cell identifiers must reference boundary data columns.');
    end
end

% Read soil type
if n_soil > 0
    if n_soil == 1
        soil(1:n_row,1:n_col)=1;                                                          % Soil type
    else
        SOIL_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["SOIL"]),2}]);
        soil = SOIL_DAT(1:n_row,1:n_col);                                                 % Soil type
        clearvars SOIL_DAT;
    end
end

% Read runoff regions
if modinf == 1 && n_runoff == 1
        runof_reg(1:n_row,1:n_col) = 1; 
elseif modinf == 1 && n_runoff > 1
        RUNOFF_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["ROFF"]),2}]);
        runof_reg = RUNOFF_DAT(1:n_row,1:n_col);                                          
        clearvars RUNOFF_DAT;
end

% Map coordinates to one cell using the raster's north-to-south row order.
x_out = floor((y_cor+n_row*s_grid-y_coord)/s_grid)+1;
y_out = floor((x_coord-x_cor)/s_grid)+1;
if n_chout > 0 && any(x_out < 1 | x_out > n_row | y_out < 1 | y_out > n_col)
    error('OFM:OutletLocation','An outlet lies outside the model grid.');
end

% Read Manning's coefficients based on Land use
n_man = str2num(InP_FILE(find(InP_FILE(:,1) == ["N_MAN"]),2));                            % The number of Manning's coefficients for roughness
if n_man == 1
    Manning = ...
        str2num(InP_FILE(find(InP_FILE(:,1) == ["MAN_VAL"]),2));                          % Manning's coefficients
    man_area(1:n_row,1:n_col) = Manning;                                                  % Manning's coefficients
else
    LAND_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["LANU"]),2}]);
    man_area = LAND_DAT(1:n_row,1:n_col);                                                 % Manning's coefficients type
    clearvars LAND_DAT;
end

% Read impervious percent
IMPE_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["IMPE"]),2}]);
imp = IMPE_DAT(1:n_row,1:n_col);                                                          % Impervious percent
clearvars IMPE_DAT;

% Read depression storage
DEPS_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["DEPS"]),2}]);
s_dep = DEPS_DAT(1:n_row,1:n_col);                                                         % Depression storage (using lumped approach, distributed approach can be developed soon)
s_dep = s_dep/100;                                                                         % Convert centimeter to meter unit.
clearvars DEPS_DAT;

active = shp ~= -9999;
roof_cells = active & (shp == 999 | dem == 999);
if any(~isfinite(dem(active))) || any(~isfinite(man_area(active)) | man_area(active) <= 0) || ...
        any(~isfinite(s_dep(active)) | s_dep(active) < 0)
    error('OFM:GridInput','Active cells need finite elevations, positive roughness and nonnegative depression storage.');
end

if modint == 1 
    % Read LAI
    LAI_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["CLAI"]),2}]);
    lai = LAI_DAT(1:n_row,1:n_col);                                                       % Read canopy leaf area index
    clearvars LAI_DAT;
end

if modets == 1 
    % Read LAI
    LAI_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["CLAI"]),2}]);
    lai = LAI_DAT(1:n_row,1:n_col);                                                      % Read canopy leaf area index
    clearvars LAI_DAT;
    % Read vegetation height
    VEGH_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["VEGH"]),2}]);
    vegheight = VEGH_DAT(1:n_row,1:n_col);                                               % Read vegetation height
    clearvars VEGH_DAT;
    % Read Canopy surface temperature
    TCAN_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["TCAN"]),2}]);
    tcan = TCAN_DAT(1:n_row,1:n_col);                                                
    clearvars TCAN_DAT;
    % Read Albedo of soil
    ALBS_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["ALBS"]),2}]);
    albs = ALBS_DAT(1:n_row,1:n_col);                                               
    clearvars ALBS_DAT;
    % Read Albedo of vegetation
    ALBV_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["ALBV"]),2}]);
    albv = ALBV_DAT(1:n_row,1:n_col);                                                   
    clearvars ALBV_DAT;
    % Read Emissivity of soil
    EMIS_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["EMIS"]),2}]);
    emis = EMIS_DAT(1:n_row,1:n_col);                                               
    clearvars EMIS_DAT;
    % Read Emissivity of vegetation
    EMIV_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["EMIV"]),2}]);
    emiv = EMIV_DAT(1:n_row,1:n_col);                                                    % Read vegetation height
    clearvars EMIV_DAT;
end

% Read initial condition
INIL_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["INIL"]),2}]);
inil = INIL_DAT(1:n_row,1:n_col);                                                         % Initial condition
clearvars INIL_DAT;

% Read rainfall (in mm)
if id_rain == 0                                                                           % Rainfall at grids
    RAIN_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["RAIN"]),2}]);
    nRain = ceil(t_dur/deltaT_rain);
    rainfall = zeros(n_row,n_col,nRain);
    for i = 1:nRain
        rainfall(:,:,i) = RAIN_DAT(1+(i-1)*n_row:i*n_row,1:n_col);
    end
else
    RAIN_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["RAIN"]),2}]);
    rainfall = RAIN_DAT(:,1:n_gauge);
end
clearvars RAIN_DAT i;

% Read weather data(forcings) for intercept model: TA (air temperature), ET (potential evaporation)
if modint == 1
    deltaT_wea = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTAT_WEA"]),2));               % Time interval for forcings (in second)
    FORC_DAT = readtable([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["FORC"]),2}],'FileType','text');
    for i=1:ceil(t_dur/deltaT_wea)
        forcingsraw(i,:) = [FORC_DAT.TA(i),FORC_DAT.RH(i), FORC_DAT.SR(i),FORC_DAT.LR(i), FORC_DAT.WD(i)];
    end  
    forcingint = forcingsraw;                         % Processing data: 0: for cummulative rainfall, 1: constant, 2: linear interpolation
    % Read Albedo of soil
    ALBS_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["ALBS"]),2}]);
    albs = ALBS_DAT(1:n_row,1:n_col);
    % Read vegetation height
    VEGH_DAT = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["VEGH"]),2}]);
    vegheight = VEGH_DAT(1:n_row,1:n_col);                                                 % Read vegetation height
    clearvars i FORC_DAT forcingsraw ALBS_DAT VEGH_DAT;
end

% Read weather data (forcings) for evapotranspiration model
if modets == 1
    deltaT_wea = str2num(InP_FILE(find(InP_FILE(:,1) == ["DELTAT_WEA"]),2));                % Time interval for forcings (in second)
    FORC_DAT = readtable([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["FORC"]),2}],'FileType','text');
    for i=1:ceil(t_dur/deltaT_wea)
        forcingsraw(i,:) = [FORC_DAT.TA(i),FORC_DAT.SR(i),FORC_DAT.WD(i),FORC_DAT.VP(i),FORC_DAT.RH(i),FORC_DAT.PA(i)];
    end 
    % Processing data: 0: for cummulative rainfall, 1: constant, 2: linear interpolation 
    forcingets(:,1) = forcingsraw(:,1);                  % Air temperature, degrees celsius
    forcingets(:,2) = forcingsraw(:,2);                  % Incoming shortwave radiation, W/m2
    forcingets(:,3) = forcingsraw(:,3);                  % Wind speed, m/s
    forcingets(:,4) = forcingsraw(:,4);                  % Air vapor pressure, kPa
    forcingets(:,5) = forcingsraw(:,5);                  % Relative humidity, kg/kg
    forcingets(:,6) = forcingsraw(:,6);                  % Atmospheric pressure, kPa
    clearvars i FORC_DAT forcingsraw;
end
              
% Parse sewer identifiers once, independent of section order and comments.
node_ad = -9999*ones(n_row,n_col);
node_invert = -9999*ones(n_row,n_col);
if moddrain == 1
    if exist('swmm_live','file') ~= 3
        error('OFM:SWMMBuild','Run build_swmm once before using MODDRA=1.');
    end
    validateattributes(DELTA_T1D,{'numeric'},{'scalar','positive','finite'});
    if deltaT < 0.001 || abs(t_dur-round(t_dur)) > 1e-8
        error('OFM:CouplingStep','SWMM requires DELTA_T >= 0.001 seconds and a whole-second run duration.');
    end
    sewerInput = [InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2}, ...
        InP_FILE{find(InP_FILE(:,1)==["DRAI"]),2}];
    [network,sewerInputText] = prepare_swmm_input(sewerInput,Start_Date,End_Date,DELTA_T1D);
    if numel(network.NodeIDs) ~= n_node || numel(network.LinkIDs) ~= n_link || ...
            numel(network.OutfallIDs) ~= n_out
        error('OFM:SWMMCount','N_NODES, N_LINKS and N_OUTS must match the SWMM network.');
    end
    nodeRows = floor((y_cor+n_row*s_grid-network.XY(:,2))/s_grid)+1;
    nodeCols = floor((network.XY(:,1)-x_cor)/s_grid)+1;
    if any(nodeRows < 1 | nodeRows > n_row | nodeCols < 1 | nodeCols > n_col)
        error('OFM:SWMMCoordinates','A coupled junction lies outside the OFM grid.');
    end
    nodeCells = sub2ind([n_row n_col],nodeRows,nodeCols);
    if numel(unique(nodeCells)) ~= n_node || any(~active(nodeCells)) || ...
            any(ismember(shp(nodeCells),[4 5 6]))
        error('OFM:SWMMCoordinates','Each junction needs its own active cell without a prescribed boundary.');
    end
    node_ad(nodeCells) = 1:n_node;
    node_invert(nodeCells) = dem(nodeCells)-network.Invert;
else
    n_node = 0; n_link = 0; n_out = 0;
end

% Read free outlet
if out_free == 1
    OUTFR = textread([InP_FILE{find(InP_FILE(:,1)==["FOLDER"]),2},InP_FILE{find(InP_FILE(:,1)==["NODF"]),2}]);
    n_ofree = OUTFR(1,1);                                                                       % The number of free outlet
    for noutf = 1:n_ofree
        x_ofree(noutf,1) = OUTFR(noutf+1,2);                                                    % Location X of free outlet)
        y_ofree(noutf,1) = OUTFR(noutf+1,3);                                                    % Location Y of free outlet
    end
    
    x_node_c = x_cor;
    y_node_c = y_cor + n_row * s_grid;
    for i=1:n_row
        for j=1:n_col
            ofree_ad(i,j) = -9999;
            x = x_node_c + j * s_grid - s_grid/2;
            y = y_node_c - i * s_grid + s_grid/2;
            for k=1:n_ofree
                if j == floor((x_ofree(k,1)-x_cor)/s_grid)+1 && i == floor((y_cor+n_row*s_grid-y_ofree(k,1))/s_grid)+1
                    ofree_ad(i,j) = k;
                end
            end
        end
    end
    clearvars i j k x y x_node_c y_node_c x_ofree y_ofree OUTFR;
end

% Load common variables
common;
display(' ');

% Report input
%
%
%
%


%% #################################################################################################################################################
%
%                                                          INITIALIZATION
%
% ##################################################################################################################################################
%% Set initial value for simulation
display(['Initializing . . .']);

vin = 0.;                                               % Total volume of rainfall
vout = 0.;
vexcess = 0.;
vinftot = 0.;                                           % Total volume of infiltration
amincdepth = 0.;
amaxcdepth = 0.;
Qpeak = zeros(1,n_chout);                                              % Maximum discharge at outlet

% For OFM
oflow = max(0,inil);
oflow(~active | shp == 1 | (shp == 2 & modchn == 1)) = 0;
%                             % Overland flow
of_rate(1:n_row,1:n_col) = 0;                           % Overland flow rate
chn_rate(1:n_row,1:n_col) = 0;                          % Channel flow rate
v_infil(1:n_row,1:n_col) = 0;                           % Infiltration depth = water table
tot_rain_depth(1:n_row,1:n_col) = 0;                    % Total rainfall depth
uvfluxs(1:n_row,1:n_col,1:2) = 0;                       % Flow rate in (u,v) directions
inun_max(1:n_row,1:n_col) = 0;                          % Initial maximum flow rate

% For RUNOFF
baseflow(1:n_row,1:n_col) = 0;                          % Baseflow
slow_interflow(1:n_row,1:n_col) = 0;                    % Slow interflow 
fast_interflow(1:n_row,1:n_col) = 0;                    % Fast interflow
sat_storage(1:n_row,1:n_col) = 0;                       % Groundwater storage
unsat_storage(1:n_row,1:n_col) = 0;                     % Upper soil storage
perc(1:n_row,1:n_col) = 0;                              % Percolation

% For CANOPY
canopystore(1:n_row,1:n_col) = 0;                       % Canopy storage
snowpackstore(1:n_row,1:n_col) = 0;                     % Snowpack storage
throughfall(1:n_row,1:n_col) = 0;                       % Throughfall
snowmetl(1:n_row,1:n_col) = 0;                          % Snowmetl

% For Infiltration
bf(1:n_row,1:n_col) = 0;                                % Cumulative precipitation for infiltration (cm)
inf_rate(1:n_row,1:n_col) = 0;                          % Infiltration rate (cm/h)
stor(1:n_row,1:n_col) = 0;                              % Soil moisture [cm]
ro(1:n_row,1:n_col) = 0;                                % Runoff (cm)
ipond(1:n_row,1:n_col) = 0;                             % Ponded conditions
fp(1:n_row,1:n_col) = 0;                                % Infiltration rate at ponding (cm/h)
tp(1:n_row,1:n_col) = 0;                                % Time to ponding (h)
tpp(1:n_row,1:n_col) = 0;                               % tpp = time required to infiltrate Fp if the system had started in ponded conditions (h)

% For Evapotranspiration
act_evap(1:n_row,1:n_col) = 0;                          % Actual Evaporation [m/s]
evap(1:n_row,1:n_col) = 0;                              % Evapotranspiration [m/s]
tot_act_evap(1:n_row,1:n_col) = 0;                      % Total evapotranspiration [m]

% For Soil water fluxes model
psim(1:n_row,1:n_col) = -500;                            % matric soil water potential, kPa
gwat(1:n_row,1:n_col) = 0;                              % groundwater storage below soil layers, mm
vrfli(1:n_row,1:n_col) = 0;                             % vertical drainage rate from layer i, mm/day
infli(1:n_row,1:n_col) = 0;                             % infiltration rate into layer, mm/day
byfli(1:n_row,1:n_col) = 0;                             % bypass flow rate from layer, mm/day
dsfli(1:n_row,1:n_col) = 0;                             % downslope flow rate from layer, mm/day
ntfli(1:n_row,1:n_col) = 0;                             % net flow rate into layer, mm/day
trani(1:n_row,1:n_col) = 0;                             % transpiration rate from layer, mm/day
swat(1:n_row,1:n_col) = 0;                              % total soil water in all layers, mm
gwfl(1:n_row,1:n_col) = 0;                              % streamflow from groundwater discharge, mm/d
seep(1:n_row,1:n_col) = 0;                              % deep seepage loss from groundwater, mm/d

% Optional soil states use a layer dimension, never a time dimension.
soil_evap=zeros(n_row,n_col); transpiration=soil_evap; potential_et=soil_evap;
soil_balance=soil_evap; seep_total=soil_evap; soil_et_total=soil_evap;
if modswf==1
    layers=modsfl_inf.nlayer;
    psim=repmat(reshape(modsfl_inf.psim_init,1,1,layers),n_row,n_col,1);
    gwat(:)=modsfl_inf.gwat_init;
    vrfli=zeros(n_row,n_col,layers); infli=vrfli; byfli=vrfli; dsfli=vrfli; ntfli=vrfli; trani=vrfli;
end

% For Drainage
node_depth_at_t = zeros(1,n_node);                          % Initial water depth of manhole
q_to_freeoutlet = [];
node_surch = []; node_depth = []; link_flow = []; out_flow = []; out_depth = [];
q_to_manhole = zeros(1,n_node);
node_flood_volume = zeros(1,n_node);
sewer_in_volume = zeros(1,n_node);
sewer_return_volume = zeros(1,n_node);
sewer_out_volume = zeros(1,n_out);

Temp_folder = fullfile(current_folder,'Results',InP_FILE{find(InP_FILE(:,1)==["RES_GEN"]),2});
% Give every run a separate directory so incomplete or older runs cannot mix.
if ~exist(Temp_folder,'dir'), mkdir(Temp_folder); end
[~,run_id] = fileparts(tempname);
Temp_folder = fullfile(Temp_folder,['run_',run_id]);
mkdir(Temp_folder);
if moddrain == 1
    drainage_folder = fullfile(Temp_folder,'drainage');
    mkdir(drainage_folder);
    preparedFile = fullfile(drainage_folder,'model.inp');
    fid = fopen(preparedFile,'wt');
    if fid < 0, error('OFM:SWMMInput','Cannot write the startup SWMM input.'); end
    fprintf(fid,'%s',sewerInputText);
    fclose(fid);
    clear sewerInputText
    sewer = drainage(preparedFile,fullfile(drainage_folder,'model.rpt'), ...
        network.NodeIDs,network.LinkIDs,network.OutfallIDs);
    sewerCleanup = onCleanup(@() sewer.close());
    node_depth_at_t = sewer.State.NodeDepth;
end
save(fullfile(Temp_folder,'metadata.mat'),'dem','shp','head_plot','deltaT', ...
    'Start_Date','End_Date','n_row','n_col','NT_RunOF','-v7.3');
display(' ');


%% #################################################################################################################################################
%
%                                                          MODEL SIMULATION
%
% ##################################################################################################################################################

% Time iterations (t)
% Snapshots are written after every computational step.
t_save = 0;
T_warning = NT_WAR;
tic;

for t = 1:NT_RunOF

    if t == 1 || t == T_warning
        display(['============ TIME STEP: ',datestr(Start_Date + seconds(deltaT) * t,'yyyy/mm/dd HH:MM:SS'),' ============']);
        display(['             Rainfall-Runoff Model is running . . .']);
    end
    
    q_to_manhole(:) = 0;
    q_to_freeoutlet(:) = 0;
    q_out = zeros(1,n_chout);
    Qout = zeros(1,n_chout);
    if n_boun > 0
        boundary_at_t = sample_forcing(boundary_dat,t,deltaT,deltaT_bound,2,false);
    end
    if modint == 1
        forcingint_at_t = sample_forcing(forcingint,t,deltaT,deltaT_wea,0,false);
    end
    if modets == 1
        forcingets_at_t = sample_forcing(forcingets,t,deltaT,deltaT_wea,0,false);
    end
    % Read rainfall
    if id_rain == 1
        RAIN_at_t = sample_forcing(rainfall,t,deltaT,deltaT_rain,1,false);                                          % Local rainfall
    elseif id_rain == 0
        RAIN_at_t = sample_forcing(rainfall,t,deltaT,deltaT_rain,1,true);                                        % Gridded rainfall
    end
    
    % Apply the rainfall to each grid cell within the watershed in order to
    % calcuate average areal rainfall at each step of rain storm
    
    % Update the boundary based on boundary condition
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) ~= -9999
                if shp(x,y) == 5                                        % Free boundary
                    oflow(x,y) = 0;
                elseif shp(x,y) == 4                                    % Soft boundary with water elevation
                    oflow(x,y) = max(0,boundary_at_t(bshape(x,y))-dem(x,y));
                end
            end
        end
    end
    % ===========================================================================================================================================
    % COMPUTE RAINFAL-RUNOFF MODEL
    % ===========================================================================================================================================  
    
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) ~= -9999
                
                % Set the value of flow rate to cell x, y to zero
                of_rate(x,y) = 0;
                
                
                % Compute rainfall depth at cell x, y
                rain_int(x,y) = rainfest(deltaT,s_grid,n_row,x,y,x_cor,y_cor,id_rain,n_gauge,x_rain,y_rain,RAIN_at_t);
                % Determine the total rainfall Depth
                tot_rain_depth(x,y) = tot_rain_depth(x,y) + rain_int(x,y)*deltaT;
                   
                % =========================================================
                % EVAPOTRANSPIRATION MODEL: Compute actual evaporation rate
                % =========================================================
                if modets == 1
                    % Use Energy balance model
                    [potential_et(x,y),evap_EB] = evapotrans(Start_Date,deltaT,t, x_cor,y_cor,n_row,x,y,s_grid,areazone,...
                                                forcingets_at_t,lai(x,y), vegheight(x,y), tcan(x,y),...
                                                albs(x,y), albv(x,y), emis(x,y), emiv(x,y),rain_int(x,y),oflow(x,y),stor(x,y),...
                                                G_daylightSavings, G_Zair, G_Zu, G_stefBoltz, G_Rv, G_cP, G_mW, G_mA, G_rhow, G_esA, G_esB, G_esC, G_kappa);
                    if modint == 0
                        evap(x,y) = evap_EB;
                    end
                end
                
                % =========================================================
                % CANOPY INTERCEPTION, SNOWMETL MODEL, and EVAPOTRANSPIRATION
                % =========================================================
                canopyPotential=[];
                if modets==1, canopyPotential=potential_et(x,y); end
                if modint == 1
                     % Use Rutter–Gash canopy formulation
                     [canopystore(x,y), snowpackstore(x,y), throughfall(x,y), snowmetl(x,y),et_soil(x,y),et_canopy(x,y),act_lum,evap(x,y),groundPotential]=intercsnow(deltaT,rain_int(x,y),...
                                                forcingint_at_t,albs(x,y),lai(x,y),vegheight(x,y),canopystore(x,y),snowpackstore(x,y),oflow(x,y),stor(x,y),...
                                                G_stefBoltz, G_kappa, G_lambdalh, G_rhoa, G_Cp, G_Rv, G_pfreez, G_Gam_ma, G_Ds, G_candec, G_F, G_zwd, G_Rc,canopyPotential);
                    % 
                    act_evap(x,y)=act_lum;
                    potential_et(x,y)=groundPotential+et_canopy(x,y);
                else
                    groundPotential=potential_et(x,y);
                    act_evap(x,y)=evap(x,y);
                end
                % Unmet ground demand can withdraw soil water only via SVAT.
                if modint==1, surfaceEvap=et_soil(x,y); else, surfaceEvap=evap(x,y); end
                unmetET=max(0,groundPotential-surfaceEvap)*(1-imp(x,y));
                
                % Compute temporary water depth from the effective rainfall and evapotranspiration
                if modint == 1
                    of_rate(x,y) = throughfall(x,y) + snowmetl(x,y) - et_soil(x,y);
                else
                    of_rate(x,y) = rain_int(x,y) - evap(x,y);
                end
                of_rate(x,y) = max(of_rate(x,y),-oflow(x,y)/deltaT);
                
                % COMPUTE INFLOW - OUTFLOW FROM POINT SOURCES -------------
                if shp(x,y) == 6                                                                % Soft boundary with water inflow [m3/s]
                    of_rate(x,y) = of_rate(x,y) + boundary_at_t(bshape(x,y)) / s_grid^2;
                end
                
                % =========================================================
                % INFILTRATION MODEL
                % =========================================================
                if modinf == 1
                    [bf(x,y),inf_rate(x,y),stor(x,y),ro(x,y),...
                        ipond(x,y),fp(x,y),tp(x,y),tpp(x,y),v_infil(x,y), of_rate(x,y)]=...
                        infiltration(bf(x,y),inf_rate(x,y),stor(x,y),ro(x,y),ipond(x,y),fp(x,y),tp(x,y),tpp(x,y),...
                        x,y,deltaT,soil,v_infil(x,y),pa_infil,oflow(x,y),of_rate(x,y),imp);
                    
                    if modswf == 0                              
                        % =========================================================
                        % RUNOFF MODEL: This MODEL runs if the SOIL WATER
                        % FLUXES MODEL is NOT selected.
                        % =========================================================
                        [of_rate(x,y),baseflow(x,y),slow_interflow(x,y),fast_interflow(x,y),sat_storage(x,y),unsat_storage(x,y),perc(x,y)] = runoffsimp...
                            (pa_runoff,runof_reg(x,y),inf_rate(x,y),sat_storage(x,y), unsat_storage(x,y),of_rate(x,y),deltaT);
                    end
                end
                
                % =========================================================
                % SOIL WATER FLUXES MODEL
                % =========================================================
                if modswf == 1
                    profile=reshape(psim(x,y,:),1,[]);
                    cover=1-exp(-0.5*max(0,lai(x,y)));
                    [profile,gwat(x,y),vFlux,iFlux,bFlux,dFlux,nFlux,tFlux,swat(x,y), ...
                        gwfl(x,y),seep(x,y),soilDiag]=soilground(modsfl_inf,profile, ...
                        cover*unmetET,gwat(x,y),(1-cover)*unmetET,inf_rate(x,y),deltaT);
                    psim(x,y,:)=reshape(profile,1,1,[]);
                    vrfli(x,y,:)=reshape(vFlux,1,1,[]); infli(x,y,:)=reshape(iFlux,1,1,[]);
                    byfli(x,y,:)=reshape(bFlux,1,1,[]); dsfli(x,y,:)=reshape(dFlux,1,1,[]);
                    ntfli(x,y,:)=reshape(nFlux,1,1,[]); trani(x,y,:)=reshape(tFlux,1,1,[]);
                    soil_evap(x,y)=soilDiag.SoilEvaporation/1000/86400;
                    transpiration(x,y)=soilDiag.Transpiration/1000/86400;
                    soil_balance(x,y)=soilDiag.BalanceResidual;
                    act_evap(x,y)=act_evap(x,y)+soil_evap(x,y)+transpiration(x,y);
                    of_rate(x,y)=of_rate(x,y)+(gwfl(x,y)+sum(bFlux+dFlux))/1000/86400;
                    seep_total(x,y)=seep_total(x,y)+seep(x,y)/1000/86400*deltaT;
                    soil_et_total(x,y)=soil_et_total(x,y)+(soil_evap(x,y)+transpiration(x,y))*deltaT;
                end
                tot_act_evap(x,y)=tot_act_evap(x,y)+act_evap(x,y)*deltaT;

                % Keeping track of the Total Volume of rainfall and total volume of infitration
                if t == NT_RunOF
                    vin = vin + tot_rain_depth(x,y)*s_grid^2;
                    vinftot = vinftot + v_infil(x,y)*s_grid^2;
                end
                
                % Compute flow to the manhole
                if moddrain == 1
                    if node_ad(x,y) ~= -9999
                        [q_to_manhole(1,node_ad(x,y)), of_rate(x,y)] = qmanhole(deltaT,s_grid,s_dep(x,y),oflow(x,y),of_rate(x,y),node_depth_at_t(node_ad(x,y)),dem(x,y),node_invert(x,y));
                    end
                end
                
                % Compute flow to the manhole
                if out_free == 1 && ofree_ad(x,y) ~= -9999
                    q_to_freeoutlet(1,ofree_ad(x,y)) = qnode(deltaT,of_rate(x,y),s_grid,s_dep(x,y),oflow(x,y));
                    of_rate(x,y) = of_rate(x,y) - q_to_freeoutlet(1,ofree_ad(x,y)) / s_grid^2;
                end
            end
        end
    end    
    
   %% ===========================================================================================================================================
    % COMPUTE OVERLAND FLOW MODEL
    % ===========================================================================================================================================
    fixedDepth = oflow;
    if TSP2D == 1
        %% ======================================================================================================================================
        % Legacy sequential routing
        % =======================================================================================================================================
        if t == 1 || t == T_warning
            display(['             Overland Flow Model is running with sequential routing . . .']);
        end
        for x = 1:n_row
            for y = 1:n_col
                if shp(x,y) ~= -9999
                    for z = -1:0                        % velocity u = (z=-1); velocity v = (z=0)
                        xx = x + z + 1;
                        yy = y - z;
                        if xx <= n_row && yy <= n_col && shp(xx,yy) ~= -9999
                            Res_flowvec = overland(TSP2D,deltaT,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,of_rate,roof_cells);
                            of_rate(x,y) = of_rate(x,y) - Res_flowvec;
                            of_rate(xx,yy) = of_rate(xx,yy) + Res_flowvec;
                        end
                    end
                end
            end
        end
        
    elseif TSP2D == 2
        %% ===========================================================================================================================================
        % Use Explicit scheme for timestepping: Adaptive timestepping solution to maintain stability of OFM
        % ===========================================================================================================================================
        if t == 1 || t == T_warning
            display(['             Overland Flow Model is running with Explicit scheme. . .']);
        end
        
        [of_rate,uvfluxs] = explicit(TSP2D,shp,n_row,n_col,deltaT,s_dep,s_grid,man_area,dem,oflow,of_rate,fixedDepth);
        
    end
    
    %%
    % ===========================================================================================================================================
    % COMPUTE DRAINAGE ROUTING
    % ===========================================================================================================================================
    drainage_updated = false;
    if moddrain == 1
        % Exchange every surface step. SWMM keeps its state and may use
        % smaller hydraulic steps; all clocks end at the same instant.
        sewerState = sewer.advanceTo(t*deltaT,q_to_manhole);
        node_flood_volume = sewerState.NodeFloodVolume;
        node_surch = node_flood_volume/deltaT;
        node_depth = sewerState.NodeDepth;
        node_depth_at_t = node_depth;
        link_flow = sewerState.LinkFlow;
        out_flow = sewerState.OutFlow;
        out_depth = sewerState.OutDepth;
        sewer_in_volume = sewer_in_volume+q_to_manhole*deltaT;
        sewer_return_volume = sewer_return_volume+node_flood_volume;
        sewer_out_volume = sewer_out_volume+sewerState.OutVolume;
        % Return this interval's actual flooded volume once, in this step.
        % Do not repeat a stale surcharge rate in the following interval.
        of_rate(nodeCells) = of_rate(nodeCells) + node_flood_volume(:)/(s_grid^2*deltaT);
        drainage_updated = true;
    end

    % ===========================================================================================================================================
    % COMPUTE CHANNEL ROUTING
    % ===========================================================================================================================================
    % If there is a routing channel => Update the channel depth
    visited_channel = false(n_row,n_col);
    if modchn == 1
        chn_rate(:) = 0;
        for i = 1:n_chn
            w_chn = pa_chn(i,1);
            d_chn = pa_chn(i,2);
            s_fact = pa_chn(i,4);
            
            for j =1:max_chn
                x = ad_chn(i,j,1);
                y = ad_chn(i,j,2);
                xx = ad_chn(i,j+1,1);
                
                % Channel depths are updated at each link and node except for the
                % last node of the link because the last node of one link is the
                % the first node of the next link
                % In the outlet link wr update depth in all nodes
                if x > 0 && y > 0 && ~visited_channel(x,y)
                    % Find new channel depth after adding inflow volume
                    vol_in = 0;
                    visited_channel(x,y) = true;
                    
                    % The flow comes from the overland (vol_ov_in)
                    if s_dep(x,y) > d_chn
                        s_dep_of = s_dep(x,y) - d_chn;
                    else
                        s_dep_of = 0;
                    end
                    vol_ov_in = 0;
                    
                    if oflow(x,y) + of_rate(x,y) * deltaT > s_dep_of
                        vol_ov_in = (oflow(x,y) + of_rate(x,y) *deltaT - s_dep_of) * s_grid^2;
                        of_rate(x,y) = (s_dep_of-oflow(x,y))/deltaT;
                        % THIS IS THE REASON WHY oflow(x,y) EQUAL TO ZERO AT OUTLET
                    end
                    tot_vol_chn = vol_in + vol_ov_in;
                    n_cdepth = cdepth(s_grid,w_chn,d_chn,s_fact,x,y,tot_vol_chn,inil);
                    inil(x,y) = n_cdepth;
                    
                    % Determine the minimum and the maximum channel depths
                    if inil(x,y) < amincdepth
                        amincdepth = inil(x,y);
                    end
                    if inil(x,y) > amaxcdepth
                        amaxcdepth = inil(x,y);
                    end
                    chn_rate(x,y) = 0;
                    if t == NT_RunOF
                        vexcess = vexcess + inil(x,y) * s_grid * w_chn;
                    end
                end
            end
        end
        
        for i = 1:n_chn
            for j =1:max_chn
                % When ad_chn(i,j+1,1) is less than zero, then that indicates that the channel
                % routing for the current link is complete. There is no routing for the outlet
                % cell
                if ad_chn(i,j,1) > 0 && ad_chn(i,j+1,1) > 0
                    [q_chn chn_rate] = chnroutin(deltaT,n_chn,s_grid,j,i,dem,s_dep,inil,ad_chn,pa_chn,chn_rate,n_flow,locat_flow);
                end
            end
        end
    end
    
    
    % Compute discharges from the surface and channel outlet storage.
    for iout = 1:n_chout
        % Calculate the flow going out from the overland portion
        pflow = sqrt(slop_out(iout)) / Man_out(iout);
        q_out(iout) = qoutlet(s_grid,deltaT,oflow,x_out,y_out,iout,of_rate,D_out,W_out,pflow);
        % Overland water depth at outlet cell is reduced after taking the outflow out of the cell
        of_rate(x_out(iout),y_out(iout)) = of_rate(x_out(iout),y_out(iout)) - q_out(iout) / s_grid^2;
        
        % Limit the channel outlet using storage after pending internal fluxes.
        if modchn == 1
            x = x_out(iout); y = y_out(iout);
            outletFactor = 1;
            for ci = 1:n_chn
                if any(ad_chn(ci,:,1) == x & ad_chn(ci,:,2) == y)
                    outletFactor = pa_chn(ci,4); break;
                end
            end
            predicted = cdepth(s_grid,W_out(iout),D_out(iout),outletFactor, ...
                x,y,chn_rate(x,y)*deltaT,inil);
            hout = max(0,predicted-s_dep(x,y));
            qchn_out = cflow(deltaT,s_grid,predicted,hout,W_out(iout), ...
                D_out(iout),s_dep(x,y),Man_out(iout),slopchn(iout),outletFactor);
            chn_rate(x,y) = chn_rate(x,y)-qchn_out;
            % Outlets absent from the channel links still need their storage updated.
            if n_chn == 0 || ~visited_channel(x,y)
                inil(x,y) = cdepth(s_grid,W_out(iout),D_out(iout),outletFactor, ...
                    x,y,chn_rate(x,y)*deltaT,inil);
                chn_rate(x,y) = 0;
            end
        else
            qchn_out = 0;
        end

        % The total outflow at the basin's outlet is given by adding the outflow from the overland & channel portion of the cell
        Qout(iout) = q_out(iout) + qchn_out;
        vout = vout + Qout(iout) * deltaT;
        if Qout(iout) > Qpeak(iout)
            Qpeak(iout) = Qout(iout);
        end
    end
    % Commit internal channel transfers in this step, once per shared node.
    if modchn == 1
        visited_channel(:) = false;
        for ci = 1:n_chn
            for cj = 1:max_chn
                x = ad_chn(ci,cj,1); y = ad_chn(ci,cj,2);
                if x > 0 && y > 0 && ~visited_channel(x,y)
                    inil(x,y) = cdepth(s_grid,pa_chn(ci,1),pa_chn(ci,2), ...
                        pa_chn(ci,4),x,y,chn_rate(x,y)*deltaT,inil);
                    visited_channel(x,y) = true;
                end
            end
        end
        chn_rate(:) = 0;
    end
    %%
    % ===========================================================================================================================================
    % UPDATE WATER DEPTH
    % ===========================================================================================================================================
    % Compute the portion of the overland depth due to: RAINFALL, INFLOW
    % AND OUTFLOW SOURCES, INFILTRATION, OVERLAND
    oflow(active) = oflow(active) + of_rate(active)*deltaT;
    oflow(shp == 4 | shp == 5) = fixedDepth(shp == 4 | shp == 5);
    if ~isreal(oflow) || any(~isfinite(oflow(active))) || any(oflow(active) < -1e-10)
        error('OFM:InvalidDepth','Invalid surface depth at step %d.',t);
    end
    oflow(active) = max(oflow(active),0);

    %% ===========================================================================================================================================
    % REPORT RESULTS
    % ===========================================================================================================================================
    if t == 1 || t == T_warning
        RunTime = toc;
        display(['             RESULT REPORTING AT TIME t = ',   num2str(t)])
        display(['             Runtime: ',   num2str(round(RunTime)),'  secs'])
        [x_ir y_ir] = find(oflow == max(oflow(:)));
        if numel(x_ir) == 1 && numel(y_ir) == 1
            
            display(['             Maximum inundation:   ',num2str(oflow(x_ir,y_ir)),'  m']);
            display(['             Grid cell:   ',num2str(x_ir),' - ',num2str(y_ir)]);
        end
        display(['             Minimum inundation:   ',num2str(min(oflow(:))),'  m']);
        display(' ');
        if t == T_warning
            T_warning = T_warning + NT_WAR;
        end
    end
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) == -9999
                inun_max(x,y) = -9999;
            elseif shp(x,y) == 2 && modchn == 1
                if inil(x,y) > inun_max(x,y)
                    inun_max(x,y) = inil(x,y);
                end
            else
                if oflow(x,y) > inun_max(x,y)
                    inun_max(x,y) = oflow(x,y);
                end
            end
        end
    end
    
    
    % Write this step only. No gridded output history is retained.
    t_save = t_save+1;
    Waterdepth = oflow;
    WaterVec = of_rate;
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) == -9999
                Results.Exp_rain(x,y) = -9999;
                Results.Exp_vol(x,y) = -9999;
                Results.Exp_flo(x,y) = -9999;
            elseif shp(x,y) == 2 && modchn == 1
                Results.Exp_rain(x,y) = rain_int(x,y);
                Results.Exp_vol(x,y) = v_infil(x,y);
                Results.Exp_flo(x,y) = inil(x,y);
            else
                Results.Exp_rain(x,y) = rain_int(x,y);
                Results.Exp_vol(x,y) = v_infil(x,y);
                Results.Exp_flo(x,y) = oflow(x,y);
            end
        end
    end
    Results.Step = t;
    Results.TimeSeconds = t*deltaT;
    Results.total_evap=tot_act_evap;
    Results.potential_et=potential_et;
    Results.DateTime = Start_Date+seconds(t*deltaT);
    Results.OutputDirectory = Temp_folder;
    Results.StepsWritten = t_save;
    Res_name = fullfile(Temp_folder,sprintf('step_%09d.mat',t));
    if modint == 1
        Results.canopystore = canopystore;
        Results.snowpackstore = snowpackstore;
        Results.throughfall = throughfall;
        Results.snowmetl = snowmetl;
        Results.act_evap(1:n_row,1:n_col) = act_evap;                   % Evapotranspiration
        Results.et_canopy(1:n_row,1:n_col) = et_canopy;                 % Canopy Evapotranspiration
        Results.et_soil(1:n_row,1:n_col) = et_soil;                     % Soil Evapotranspiration
    end
    if modets == 1
        Results.act_evap(1:n_row,1:n_col) = act_evap;                   % Evapotranspiration
    end
    if modinf == 1
        Results.ro = ro;                                                % Runoff potential
        Results.stor = stor;                                            % Soil moisture
        Results.inf_rate = inf_rate;                                    % Infiltration rate
        Results.v_infil = v_infil;                                      % Water table depth
        Results.baseflow = baseflow;                                    % Baseflow
        Results.slow_interflow = slow_interflow;                        % Slow interflow
        Results.fast_interflow = fast_interflow;                        % Fast interflow
        Results.sat_storage = sat_storage;                              % Groundwater storage
        Results.unsat_storage = unsat_storage;                          % Upper soil storage
        Results.perc = perc;                                            % Percolation
    end
    if moddrain == 1
        Results.q_to_manhole = q_to_manhole;
        Results.node_surch = node_surch;
        Results.node_depth = node_depth;
        Results.link_flow = link_flow;
        Results.out_flow = out_flow;
        Results.out_depth = out_depth;
        Results.drainage_updated = drainage_updated;
        Results.SewerTimeSeconds = sewerState.TimeSeconds;
        Results.NodeFloodVolume = node_flood_volume;
        Results.SewerInVolume = sewer_in_volume;
        Results.SewerReturnVolume = sewer_return_volume;
        Results.SewerOutVolume = sewer_out_volume;
        Results.SewerVolume = sewerState.SewerVolume;
    end
    if modswf == 1
        Results.soil_evap=soil_evap;
        Results.transpiration=transpiration;
        Results.soil_balance=soil_balance;
        Results.seep_total=seep_total;
        Results.soil_et_total=soil_et_total;
        Results.psim = psim;                                            % matric soil water potential, kPa
        Results.gwat = gwat;                                            % groundwater storage below soil layers, mm
        Results.vrfli = vrfli;                                          % vertical drainage rate from layer i, mm/day
        Results.infli = infli;                                          % infiltration rate into layer, mm/day
        Results.byfli = byfli;                                          % bypass flow rate from layer, mm/day
        Results.dsfli = dsfli;                                          % downslope flow rate from layer, mm/day
        Results.ntfli = ntfli;                                          % net flow rate into layer, mm/day
        Results.trani = trani;                                          % transpiration rate from layer, mm/day
        Results.swat = swat;                                            % total soil water in all layers, mm
        Results.gwfl = gwfl;                                            % streamflow from groundwater discharge, mm/d
        Results.seep = seep;                                            % deep seepage loss from groundwater, mm/d
    end
    if modchn == 1, Results.ChannelDepth = inil; end
    Results.q_to_freeoutlet = q_to_freeoutlet;
    Results.q_out = q_out;
    Results.Qout = Qout;
    save(Res_name,'Results','Waterdepth','WaterVec','-v7.3');
end

if moddrain == 1, sewer.close(); end
FinalTime = toc;
Res_name = fullfile(Temp_folder,'summary.mat');
save(Res_name,'Results','inun_max','dem','oflow','head_plot','deltaT', ...
    'n_row','n_col','FinalTime','Qpeak','vin','vout','vinftot','-v7.3');
save(fullfile(Temp_folder,'MaxInundation.txt'),'inun_max','-ascii');
if options.Plot
    visualize(inun_max,head_plot,t,deltaT,dem,n_row,n_col);
end
% Restore the caller's directory even if the simulation raises an error.

end

%% #################################################################################################################################################
%
%                                                          SIMULATION ENDDED
%
% ##################################################################################################################################################
