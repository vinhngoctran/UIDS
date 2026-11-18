function run_swmm_tests
%RUN_SWMM_TESTS Exercise the real native engine and full OFM/SWMM coupling.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
if exist('swmm_live','file') ~= 3, error('Run build_swmm before run_swmm_tests.'); end
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
source = fullfile(folder,'network.inp');
text = strjoin(["[TITLE]"; "OFM direct exchange regression"; "[OPTIONS]"; ...
    "FLOW_UNITS CMS"; "FLOW_ROUTING DYNWAVE"; "START_DATE 01/01/2016"; ...
    "START_TIME 00:00:00"; "END_DATE 01/01/2016"; "END_TIME 00:01:00"; ...
    "REPORT_STEP 00:00:01"; "WET_STEP 00:00:01"; "DRY_STEP 00:00:01"; ...
    "ROUTING_STEP 0.1"; "VARIABLE_STEP 0.75"; "MINIMUM_STEP 0.001"; ...
    "ALLOW_PONDING YES"; "MIN_SURFAREA 0.001"; "[JUNCTIONS]"; ...
    "J2 -0.1 0.2 0.01 0 10"; "J1 0 0.05 0.02 0 10"; ...
    "[OUTFALLS]"; "O1 -0.2 FREE NO"; "[CONDUITS]"; ...
    "C1 J1 J2 5 0.013 0 0 0"; "C2 J2 O1 5 0.013 0 0 0"; ...
    "[XSECTIONS]"; "C1 CIRCULAR 0.03 0 0 0 1"; "C2 CIRCULAR 0.03 0 0 0 1"; ...
    "[COORDINATES]"; "J1 1 3"; "O1 7 3"; "J2 5 3"; ...
    "[INFLOWS]"; "J1 FLOW OLD FLOW 1 1 0"; ...
    "[TIMESERIES]"; "OLD 0:00 100"; "OLD 0:01 100"; ...
    "[FILES]"; "SAVE HOTSTART obsolete.hot"],newline);
write_text(source,text);
startTime = datetime(2016,1,1);
[network,prepared] = prepare_swmm_input(source,startTime,startTime+seconds(20),0.01);
assert(isequal(network.NodeIDs,{'J2','J1'}));
assert(isequal(network.XY,[5 3; 1 3]));
assert(~contains(prepared,'HOTSTART') && ~contains(prepared,'J1 FLOW OLD'));
input = fullfile(folder,'prepared.inp'); write_text(input,prepared);
report = fullfile(folder,'model.rpt');
% Deliberately reverse identifiers to prove values map by ID, not file order.
session = drainage(input,report,{'J1','J2'},{'C2','C1'},{'O1'});
sessionCleanup = onCleanup(@() session.close());
assert(session.State.TimeSeconds == 0 && session.State.EngineVersion == 52004);
assert(max(abs(session.State.NodeDepth-[0.02 0.01])) < 1e-12);
initial = session.State.SewerVolume;
% File sizes and mtimes must stay unchanged during the exchange loop.
before = file_state(folder);
flood = 0; out = 0; injected = 0;
for dt = [0.1 0.15 0.25 0.5 1 2]
    state = session.advance(dt,[0.05 0]);
    injected = injected+0.05*dt;
    flood = flood+sum(state.NodeFloodVolume);
    out = out+sum(state.OutVolume);
    assert(abs(state.TimeSeconds-injected/0.05) < 1e-9);
    assert(all(isfinite(state.NodeDepth)) && all(state.NodeFloodVolume >= 0));
end
assert(isequal(before,file_state(folder)),'A synchronization step wrote a file.');
assert(flood > 0,'Fixture must exercise sewer-to-surface flooding.');
fprintf('SWMM volume balance: stored %.12g, flood %.12g, outlet %.12g, injected %.12g [m^3]\n', ...
    state.SewerVolume-initial,flood,out,injected);
assert(abs(state.SewerVolume-initial+flood+out-injected) < 0.02*injected);
% Turning inflow off must not replay the preceding interval's pulse.
state = session.advance(0.1,[0 0]);
assert(abs(state.TimeSeconds-4.1) < 1e-9);
session.close(); clear sessionCleanup
% Reopen starts a new simulation, and close is idempotent.
session = drainage(input,report,{'J1','J2'},{'C1','C2'},{'O1'});
assert(session.State.TimeSeconds == 0);
session.close(); session.close();
% A missing ID must release the failed session, allowing the next run.
try
    bad = drainage(input,report,{'missing'},{'C1'},{'O1'}); %#ok<NASGU>
    error('OFM:Test','Missing ID was accepted.');
catch err
    assert(strcmp(err.identifier,'OFM:SWMMObject'));
end
session = drainage(input,report,{'J1','J2'},{'C1','C2'},{'O1'});
session.close();

% Full coupled run: surface+sewer exchange closes at every surface step.
base = fileread(fullfile(root,'Test_01.inf'));
for name = ["DEM.txt","MASK.txt","DEPRESION.txt","IMPERVIOUS.txt"]
    writematrix(zeros(3,4),fullfile(folder,name),'Delimiter','tab');
end
write_text(fullfile(folder,'HEADING.txt'),sprintf( ...
    'UTM 52 N\nncols 4\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 2\nnodata -9999\n'));
writematrix([100;100],fullfile(folder,'RAINFALL.txt'),'Delimiter','tab');
values = {'FOLDER',[folder filesep]; 'START_TIME','2016 1 1 0 0 0'; ...
    'END_TIME','2016 1 1 0 0 2'; 'DELTA_T','0.1'; 'DELTA_T1D','2'; ...
    'DELTAT_RAIN','1'; 'N_BOUN','0'; 'MODDRA','1'; 'N_NODES','2'; ...
    'N_LINKS','2'; 'N_OUTS','1'; 'DRAI','network.inp'; ...
    'RES_GEN',['regression_' char(java.util.UUID.randomUUID)]};
for k = 1:size(values,1)
    base = regexprep(base,['(?m)^' values{k,1} '[ \t]+[^\r\n]*'], ...
        [values{k,1} ' ' values{k,2}]);
end
config = fullfile(folder,'model.inf'); write_text(config,base);
old = pwd; restore = onCleanup(@() cd(old)); cd(root);
result = ofm_urban(config);
outputCleanup = onCleanup(@() rmdir(fileparts(result.OutputDirectory),'s'));
assert(result.StepsWritten == 20 && result.SewerTimeSeconds == 2);
assert(sum(result.SewerInVolume) > 0 && sum(result.SewerReturnVolume) > 0);
for step = 1:20
    saved = load(fullfile(result.OutputDirectory,sprintf('step_%09d.mat',step)),'Results');
    r = saved.Results;
    expectedRain = step*0.1*0.1*12*4;
    balance = sum(r.Exp_flo,'all')*4+sum(r.SewerInVolume)-sum(r.SewerReturnVolume);
    assert(abs(balance-expectedRain) < 1e-10);
    assert(abs(r.SewerTimeSeconds-r.TimeSeconds) < 1e-9);
    assert(r.drainage_updated);
end
assert(isempty(dir(fullfile(result.OutputDirectory,'drainage','*HS*'))));
assert(isempty(dir(fullfile(result.OutputDirectory,'drainage','*.out'))));
fprintf('All direct SWMM and coupled exchange checks passed.\n');
end

function write_text(filename,text)
fid = fopen(filename,'w'); assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function state = file_state(folder)
files = dir(folder); files = files(~[files.isdir]);
state = [{files.name};num2cell([files.bytes]);num2cell([files.datenum])];
end
