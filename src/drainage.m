%============================================================================
%   drainage.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is to simulate drainage swmm model for sewer system
%============================================================================

function [node_surch, node_depth, link_flow, out_flow, out_depth] = drainage(n_node, n_out,DRAI_DAT,Start_Date,DELTA_T1D,deltaT,t,t1D,q_to_manhole_t)

%% MODEL VARIABLES

% SWMM Type
swmm.Subcatchments = 0;
swmm.Nodes         = 1;
swmm.Links         = 2;
swmm.System        = 3;

% SWMM Units
swmm.CFS           = 0;
swmm.GPM           = 1;
swmm.MGD           = 2;
swmm.CMS           = 3;
swmm.LPS           = 4;
swmm.LPD           = 5;

%  subcatchment variables
swmm.subcatch.rainfall     = 0;  %rainfall (in/hr or mm/hr)
swmm.subcatch.snowdepth    = 1;  %snow depth (in or mm)
swmm.subcatch.loss         = 2;  %evaporation loss (in/day or mm/day)
swmm.subcatch.infiltration = 3;  %infiltration losses (in/hr or mm/hr)
swmm.subcatch.runoff       = 4;  %runoff rate (flow units)
swmm.subcatch.gwoutflow    = 5;  %groundwater outflow rate (flow units)
swmm.subcatch.gwelevation  = 6;  %groundwater water table elevation (ft or m)
swmm.subcatch.unsaturated  = 7;  %unsaturated zone moisture content (fraction)
swmm.subcatch.con1stpoll   = 8;  %runoff concentration of first pollutant

%  node variables
swmm.node.depth           = 0;  %depth of water above invert (ft or m)
swmm.node.head            = 1;  %hydraulic head (ft or m)
swmm.node.volume          = 2;  %volume of stored + ponded water (ft3 or m3)
swmm.node.lateralinflow   = 3;  %lateral inflow (flow units)
swmm.node.totalinflow     = 4;  %total inflow (lateral + upstream) (flow units)
swmm.node.flooding        = 5;  %flow lost to flooding (flow units)
swmm.node.con1stpoll      = 6;  %concentration of first pollutant

%  link variables
swmm.link.flow            = 0;  %flow rate (flow units)
swmm.link.depth           = 1;  %flow depth (ft or m)
swmm.link.velocity        = 2;  %flow velocity (ft/s or m/s)
swmm.link.volume          = 3;  %flow volume (ft3 or m3)
swmm.link.filledfrac      = 4;  %fraction of conduit's area filled or settingnon-conduits
swmm.link.con1stpoll      = 5;  %concentration of first pollutant

%  system-wide variables
swmm.sys.temp             = 0;  %air temperature (deg. F or deg. C)
swmm.sys.rainfall         = 0;  %rainfall (in/hr or mm/hr)
swmm.sys.snowdepth        = 0;  %snow depth (in or mm)
swmm.sys.lossrate         = 0;  %evaporation + infiltration loss rate (in/hr or mm/hr)
swmm.sys.runoff           = 0;  %runoff flow (flow units)
swmm.sys.dryweatherinflow = 0;  %dry weather inflow (flow units)
swmm.sys.gwinflow         = 0;  %groundwater inflow (flow units)
swmm.sys.RDIIinflow       = 0;  %RDII inflow (flow units)
swmm.sys.directinflow     = 0;  %user supplied direct inflow (flow units)
swmm.sys.totallatinflow   = 0;  %total lateral inflow (sum of variables 4 to 8) (flow units)
swmm.sys.flooding         = 0;  %flow lost to flooding (flow units)
swmm.sys.outfalls         = 0;  %flow leaving through outfalls (flow units)
swmm.sys.storedvolume     = 0;  %volume of stored water (ft3 or m3)
swmm.sys.actevap          = 0;  %actual evaporation rate (in/day or mm/day)
swmm.sys.potevap          = 0;  %potential evaporation rate (PET) (in/day or mm/day)

%% Re-organize inflow to manhole according to DELTA_T1D
if DELTA_T1D == deltaT
    q_to_manhole = q_to_manhole_t;
else
    for dt = 1:t1D
        q_to_manhole(dt,:) = mean(q_to_manhole_t((dt-1)*round(DELTA_T1D/deltaT)+1:dt*round(DELTA_T1D/deltaT),:));
    end
end

% Update input for Drainage model
% Update simulation time
NEW_DRAI_DAT = DRAI_DAT;
idx1 = find(NEW_DRAI_DAT(:,1) == ["[TIMESERIES]"]);
idx2 = find(NEW_DRAI_DAT(:,1) == ["[INFLOWS]"]);
if t1D == 1
        
    idx = find(NEW_DRAI_DAT(:,1) == ["END_DATE"]);
    NEW_DRAI_DAT{idx,2} =   datestr(Start_Date + seconds(deltaT) * t,'mm/dd/yyyy');               % End Date   
    idx = find(NEW_DRAI_DAT(:,1) == ["END_TIME"]);
    NEW_DRAI_DAT{idx,2} =   datestr(Start_Date + seconds(deltaT) * t,'HH:MM:SS');                 % End Time
    idx = find(NEW_DRAI_DAT(:,1) == '[FILES]');
    NEW_DRAI_DAT{idx+2,1} = 'SAVE';
    NEW_DRAI_DAT = fillmissing(NEW_DRAI_DAT,'constant',"");
    
    % Write input file
    PCUW1D = fopen('PCUW1D.mdf','wt');
    for iii = 1:idx1+2
        fprintf(PCUW1D,'%s\t',NEW_DRAI_DAT(iii,:));
        fprintf(PCUW1D,'\n');
    end
    % Write inflow to manhole: TimeSerie Name / Time / Inflow [cms]
    for nod = 1: n_node
        fprintf(PCUW1D,'%s\t',NEW_DRAI_DAT(idx2+2+nod,3),datestr(Start_Date + seconds(deltaT) * t,'mm-dd-yyyy HH:MM'),num2str(round(q_to_manhole(t1D,nod),5)));
        fprintf(PCUW1D,'\n');
    end
    fclose(PCUW1D);
    
else
    OLD_DATE = ['END_DATE	',datestr(Start_Date + seconds(DELTA_T1D) * (t1D-1),'mm/dd/yyyy')];         % Old End Time
    NEW_DATE = ['END_DATE	',datestr(Start_Date + seconds(deltaT) * t,'mm/dd/yyyy')];                  % New End Time 
    OLD_Time = ['END_TIME	',datestr(Start_Date + seconds(DELTA_T1D) * (t1D-1),'HH:MM:SS')];           % Old End Time
    NEW_Time = ['END_TIME	',datestr(Start_Date + seconds(deltaT) * t,'HH:MM:SS')];                    % New End Time
    Hotstart = 'USE/SAVE	HOTSTART	PCUW1D_HS';
    PCUW1D = fileread('PCUW1D.mdf');
    PCUW1D = regexprep(PCUW1D, OLD_DATE, sprintf(NEW_DATE));
    PCUW1D = regexprep(PCUW1D, OLD_Time, sprintf(NEW_Time));
    if t1D==2
        PCUW1D = regexprep(PCUW1D, 'SAVE	HOTSTART	PCUW1D_HS', sprintf(Hotstart));                 % Update Hotstart
    end
    fid = fopen('PCUW1D.mdf', 'w');
    fwrite(fid, PCUW1D);
    fclose(fid);
    
    % Write inflow to manhole: TimeSerie Name / Time / Inflow [cms]
    fid = fopen('PCUW1D.mdf','a+');
    for nod = 1: n_node
        fprintf(fid, '\n');
        fprintf(fid,'%s\t',NEW_DRAI_DAT(idx2+2+nod,3),datestr(Start_Date + seconds(deltaT) * t,'mm-dd-yyyy HH:MM'),num2str(round(q_to_manhole(t1D,nod),5)));
    end
    fclose(fid);
end

%% MODEL SIMULATION
infile = 'PCUW1D.mdf';
if ~(libisloaded('swmm5'))
    loadlibrary('swmm5');
end
rpt_file = strrep(lower(infile), '.mdf', '.rpt');
out_file = strrep(lower(infile), '.mdf', '.out');
error = calllib('swmm5','swmm_open',infile, rpt_file, out_file);
if ~ismember(1, [0, 1])
    if libisloaded('swmm5')
        unloadlibrary swmm5;
    end
end
if ~(libisloaded('swmm5'))
    loadlibrary('swmm5');
end
error = calllib('swmm5','swmm_start', 1);
elapsed_time = 1e-6;
timePtr = libpointer('doublePtr', elapsed_time);
while ~timePtr.value == 0
    if ~(libisloaded('swmm5'))
        loadlibrary('swmm5');
    end
    error = calllib('swmm5','swmm_step', timePtr);
    time  = timePtr.value*24;
end
if ~(libisloaded('swmm5'))
    loadlibrary('swmm5');
end
error = calllib('swmm5','swmm_end');
if ~(libisloaded('swmm5'))
    loadlibrary('swmm5');
end
error = calllib('swmm5','swmm_report');
if ~(libisloaded('swmm5'))
    loadlibrary('swmm5');
end
% error = calllib('swmm5','swmm_save_results');         % Save to folders
error = calllib('swmm5','swmm_close');
if libisloaded('swmm5')
    unloadlibrary swmm5;
end

%% RESULTS
% Read information results
infor = read1D(out_file);

% Surcharge
pnode_surch = readresults(out_file,infor,swmm.Nodes,swmm.node.flooding);
node_surch = pnode_surch(end,1:n_node);

% Water depth at manhole
pnode_depth = readresults(out_file,infor,swmm.Nodes,swmm.node.depth);
node_depth = pnode_depth(end,1:n_node);

% Flow in the links
plink_flow = readresults(out_file,infor,swmm.Links,swmm.link.flow);
link_flow = plink_flow(end,:);

% Volume at the outlet
pout_flow = readresults(out_file,infor,swmm.Nodes,swmm.node.totalinflow);
out_flow = pout_flow(end,end-n_out+1:end);

% Water depth at outlet
out_depth = pnode_depth(end,end-n_out+1:end);

end

%% Get time series results
function value = readresults(outFile,parameters,outtype,outid )
SUBCATCH   = 0;
NODE       = 1;
LINK       = 2;
SYS        = 3;
RECORDSIZE = 4;

N_period        = parameters(1);
FlowUnits       = parameters(2);
N_bas           = parameters(3);
N_nod           = parameters(4);
N_link          = parameters(5);
N_pollut        = parameters(6);
StartDate       = parameters(7);
ReportStep      = parameters(8);
SubcatchVars    = parameters(9);
NodeVars        = parameters(10);
LinkVars        = parameters(11);
SysVars         = parameters(12);
StartPos        = parameters(13);
BytesPerPeriod  = parameters(14);

switch outtype
    case 0
        Nval=N_bas;
    case 1
        Nval=N_nod;
    case 2
        Nval=N_link;
    otherwise
        Nval=-999;
end
if ~( Nval==-999)
    Fout = fopen(outFile, 'r');
    value = zeros(N_period,Nval);
    for j=1:Nval
        for i=1:N_period        
            offset = StartPos + (i-1)*BytesPerPeriod + 2*RECORDSIZE;
            if ( outtype == SUBCATCH )
                offset = offset + RECORDSIZE*((j-1)*SubcatchVars + outid);
            elseif (outtype == NODE)
                offset = offset + RECORDSIZE*(N_bas*SubcatchVars + (j-1)*NodeVars + outid);
            elseif (outtype == LINK)
                offset = offset + RECORDSIZE*(N_bas*SubcatchVars + N_nod*NodeVars + (j-1)*LinkVars + outid);
            elseif (outtype == SYS)
                offset = offset + RECORDSIZE*(N_bas*SubcatchVars + N_nod*NodeVars + N_link*LinkVars + outid);
            else
            end          
            SEEK_SET = fseek(Fout, offset, 'bof');
            value(i,j) = fread(Fout,1,'real*4');
        end
    end
    fclose(Fout);
end
end

