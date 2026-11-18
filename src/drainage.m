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

%% Run one interval. Only its mean surface inflow is retained by the caller.
if exist('read1D','file') ~= 3 && exist('read1D','file') ~= 2
    error('OFM:SWMMDependency','A read1D binary for this MATLAB platform is required. The repository supplies Windows binaries.');
end
intervalStart = Start_Date+seconds(t*deltaT-DELTA_T1D);
if size(q_to_manhole_t,1) == 1
    inflow = q_to_manhole_t;
else
    count = round(DELTA_T1D/deltaT);
    inflow = mean(q_to_manhole_t(end-count+1:end,:),1);
end
infile = 'PCUW1D.mdf';
rpt_file = 'PCUW1D.rpt';
out_file = 'PCUW1D.out';
write_swmm_interval(infile,DRAI_DAT,n_node,inflow,intervalStart,DELTA_T1D,t1D);
if ~libisloaded('swmm5'), loadlibrary('swmm5'); end
cleanup = onCleanup(@close_swmm);
check_swmm(calllib('swmm5','swmm_open',infile,rpt_file,out_file),'open');
check_swmm(calllib('swmm5','swmm_start',1),'start');
timePtr = libpointer('doublePtr',1e-6);
while timePtr.Value > 0
    check_swmm(calllib('swmm5','swmm_step',timePtr),'step');
end
check_swmm(calllib('swmm5','swmm_end'),'end');
check_swmm(calllib('swmm5','swmm_report'),'report');
clear cleanup

%% RESULTS
% Read information results
infor = read1D(out_file);

% Surcharge
pnode_surch = readresults(out_file,infor,swmm.Nodes,swmm.node.flooding,true);
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
function value = readresults(outFile,parameters,outtype,outid,average)
if nargin < 5, average = false; end
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
    if Fout < 0, error('OFM:SWMMOutput','Cannot open %s.',outFile); end
    cleanup = onCleanup(@() fclose(Fout));
    if N_period < 1, error('OFM:SWMMOutput','SWMM produced no report records.'); end
    value = zeros(1,Nval);
    periods = N_period;
    if average, periods = 1:N_period; end
    for j=1:Nval
        for i=periods
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
            value(1,j) = value(1,j)+fread(Fout,1,'real*4')/numel(periods);
        end
    end
    clear cleanup
end
end


function check_swmm(code,operation)
if code ~= 0
    error('OFM:SWMM','SWMM %s failed with error code %d.',operation,code);
end
end

function close_swmm
if libisloaded('swmm5')
    calllib('swmm5','swmm_close');
    unloadlibrary('swmm5');
end
end
