function run_regression_tests
%RUN_REGRESSION_TESTS Numerical invariants and short model/output regressions.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
old = pwd; cleanup = onCleanup(@() cd(old)); cd(root);

% Closed-domain volume, symmetry, sources, depression storage and no-data.
shp = zeros(5); z = zeros(5); n = 0.04*ones(5); dep = zeros(5);
h = zeros(5); h(3,3) = 0.3;
[r,f] = explicit2(2,shp,5,5,2,dep,2,n,z,h,zeros(5));
next = h+2*r;
assert(all(next(:) >= 0));
assert(abs(sum(next(:))-sum(h(:))) < 1e-12);
assert(max(abs(next-flipud(next)),[],'all') < 1e-12);
assert(max(abs(next-fliplr(next)),[],'all') < 1e-12);
assert(isequal(size(f),[5 5 2]));
[r,~] = explicit2(2,shp,5,5,2,dep,2,n,z,zeros(5),ones(5)*1e-4);
assert(max(abs(r(:)-1e-4)) < 1e-12);
dep(:) = 0.05; h(:) = 0; h(3,3) = 0.3;
[r,~] = explicit2(2,shp,5,5,20,dep,2,n,z,h,zeros(5));
next = h+20*r;
assert(next(3,3) >= 0.05-1e-12);
assert(abs(sum(next(:))-0.3) < 1e-12);
z = [0 2]; h = [0.5 0]; dep = [0 0]; n = [0.03 0.03];
assert(overland2(2,1,dep,1,n,z,h,1,1,1,2) == 0);
z = [1200 1200.2]; h = [0.4 0.2];
assert(abs(overland2(2,1,dep,1,n,z,h,1,1,1,2)) < 1e-10);
h = [0.4 0.1];
q1 = overland2(2,1,dep,1,n,z,h,1,1,1,2);
q2 = overland2(2,1,dep,1,n,z,h,1,2,1,1);
assert(abs(q1+q2) < 1e-12 && q1 > 0);
shp = [0 -9999]; n(2) = -9999;
[r,~] = explicit2(2,shp,1,2,1,dep,1,n,z,h,[0 0]);
assert(all(r == 0));

% Native forcing preserves cumulative totals with non-integer step ratios.
x = processing([6;12;3],4,6,12,1);
assert(max(abs(x-[4;6;8])) < 1e-12);
x = processing([0;6;12],4,6,12,2);
assert(max(abs(x-[4;8;12])) < 1e-12);
grid = ones(2,3)*6;
assert(isequal(sample_forcing(grid,1,2,6,1,true),ones(2,3)*2));
assert(rainfest(1,1,1,1,1,0,0,1,2,[0;1],[0.5;0.5],[1 2]) == 0.001);

% Infiltration handles dry, impermeable and ponded conditions without NaNs.
p = [1 20 0.45 0.2 0.5];
[F,rate,~,~,~,~,~,~,v,r] = infiltration(0,0,0,0,0,0,0,0,1,1,60,1,0,p,0,0,0);
assert(F == 0 && rate == 0 && v == 0 && r == 0);
[F,rate,~,~,~,~,~,~,v,r] = infiltration(0,0,0,0,0,0,0,0,1,1,60,1,0,p,0,1e-6,0);
assert(abs(rate-1e-6) < 1e-12 && abs(v-60e-6) < 1e-12);
assert(abs(F/100-v) < 1e-12 && abs(r) < 1e-12);
[F,rate,~,~,~,~,~,~,v,r] = infiltration(0,0,0,0,0,0,0,0,1,1,60,1,0,p,0.1,0,0.5);
assert(rate > 0 && v <= 0.05 && abs(v+60*r) < 1e-12);
% Integrated ponded Green-Ampt equation, F is on the pervious fraction [cm].
A = p(2)*(p(3)-p(4));
assert(abs(F-A*log1p(F/A)-p(1)/60) < 1e-10);
[~,rate,~,~,~,~,~,~,v,~] = infiltration(0,0,0,0,0,0,0,0,1,1,60,1,0,p,1,0,1);
assert(rate == 0 && v == 0);

% Reservoir conservation across short and long steps, including full drainage.
for dt = [0.5 60 8640000]
    [r,b,s,f,rs,ls,pc] = runoffsimp([1 2 3 4 0.1 0.03],1,1e-5,0.1,0.2,0,dt);
    assert(abs(rs+ls+dt*(b+s+f)-(0.3+dt*1e-5)) < 1e-10);
    assert(all([r b s f rs ls pc] >= 0));
end
for dt = [0.1 1 1.5 60]
    q = qnode(dt,0,2,0.01,0.1);
    assert(q >= 0 && q*dt <= 0.09*4+1e-12);
    [q,r] = qmanhole(dt,2,0.01,0.1,0,0,10,2);
    assert(q >= 0 && 0.1+r*dt >= 0.01-1e-12);
    q = qoutlet(2,dt,0.1,1,1,1,0,0.01,2,10);
    assert(q >= 0 && q*dt <= 0.09*4+1e-12);
end
assert(abs(cdepth(10,2,1,1.2,1,1,24,0)-1) < 1e-12);
assert(abs(cdepth(10,2,1,1.2,1,1,144,0)-2) < 1e-12);
q = cflow(100,2,1,1,1,2,0,0.03,-1);
assert(q < 0 && abs(q)*100 <= 2+1e-12);

% End-to-end run with real example input layout, including disk snapshots.
fixture = tempname; mkdir(fixture);
fixtureCleanup = onCleanup(@() rmdir(fixture,'s'));
base = fileread(fullfile(root,'Test_01.inf'));
write_matrix(fullfile(fixture,'DEM.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'MASK.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'DEPRESION.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'IMPERVIOUS.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'RAINFALL.txt'),ones(3,1));
fid = fopen(fullfile(fixture,'HEADING.txt'),'w');
fprintf(fid,'UTM 52 N\nncols 4\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 2\nnodata -9999\n'); fclose(fid);
base = setkey(base,'FOLDER',[fixture filesep]);
base = setkey(base,'END_TIME','2016 7 5 0 0 6');
base = setkey(base,'DELTA_T','0.5');
base = setkey(base,'DELTAT_RAIN','2');
base = setkey(base,'N_BOUN','0');
base = setkey(base,'RES_GEN',['regression_' char(java.util.UUID.randomUUID)]);
config = fullfile(fixture,'model.inf');
for scheme = [1 2]
    text = setkey(base,'TSTP2D',num2str(scheme));
    fid = fopen(config,'w'); fprintf(fid,'%s',text); fclose(fid);
    result = ofm_urban(config);
    outputCleanup = onCleanup(@() rmdir(fileparts(result.OutputDirectory),'s'));
    files = dir(fullfile(result.OutputDirectory,'step_*.mat'));
    assert(numel(files) == 12 && result.StepsWritten == 12);
    assert(isequal(size(result.Exp_flo),[3 4]));
    assert(max(abs(result.Exp_flo(:)-0.003)) < 1e-12);
    first = load(fullfile(result.OutputDirectory,'step_000000001.mat'));
    last = load(fullfile(result.OutputDirectory,'step_000000012.mat'));
    assert(isequal(size(first.Results.Exp_flo),size(last.Results.Exp_flo)));
    assert(isequal(size(last.Waterdepth),[3 4]));
    assert(last.Results.TimeSeconds == 6 && first.Results.TimeSeconds == 0.5);
    [time,point] = read_ofm_series(result.OutputDirectory,[1 1; 3 4],[1 12]);
    assert(isequal(time,[0.5;6]) && max(abs(point(2,:)-0.003)) < 1e-12);
    summary = load(fullfile(result.OutputDirectory,'summary.mat'));
    assert(max(abs(summary.inun_max(:)-0.003)) < 1e-12);
    clear outputCleanup
end
% Initial conditions and infiltration/runoff coupling retain all water.
write_matrix(fullfile(fixture,'DEPRESION.txt'),ones(3,4)*0.02);
text = setkey(base,'TSTP2D','2');
text = setkey(text,'MODINF','1');
result = run_case(text,fixture);
stored = result.Exp_flo+result.sat_storage+result.unsat_storage;
assert(max(abs(stored(:)-0.023)) < 1e-11);
assert(all(result.v_infil(:) > 0));

% Prescribed stages and free boundaries remain fixed at the step endpoint.
mask = zeros(3,4); mask(:,1) = 4; mask(:,4) = 5;
write_matrix(fullfile(fixture,'MASK.txt'),mask);
write_matrix(fullfile(fixture,'BOUNDARY_SHAPE.txt'),ones(3,4));
write_matrix(fullfile(fixture,'BOUNDARY.txt'),[0.02;0.03;0.04;0.05]);
text = setkey(base,'N_BOUN','1');
text = setkey(text,'DELTAT_BOUN','2');
for scheme = [1 2]
    text = setkey(text,'TSTP2D',num2str(scheme));
    result = run_case(text,fixture);
    assert(max(abs(result.Exp_flo(:,1)-0.05)) < 1e-12);
    assert(all(result.Exp_flo(:,4) == 0));
end

% Canopy with zero LAI avoids 0/0 and closes the rain/snow/storage balance.
write_matrix(fullfile(fixture,'MASK.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'DEPRESION.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'TestLAI.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'VegHeight.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'albsoil.txt'),ones(3,4)*0.2);
fid = fopen(fullfile(fixture,'weather.mdf'),'w');
fprintf(fid,'TA RH SR LR WD\n20 0.005 200 300 1\n'); fclose(fid);
text = setkey(base,'MODINT','1');
result = run_case(text,fixture);
assert(all(isfinite(result.Exp_flo),'all'));
assert(all(result.canopystore == 0,'all'));
assert(all(result.et_canopy == 0,'all'));

% Channel terminal cells and reverse transfers conserve storage this step.
mask = zeros(3,4); mask(2,1:3) = 1;
initial = zeros(3,4); initial(2,1) = 0.2;
write_matrix(fullfile(fixture,'MASK.txt'),mask);
write_matrix(fullfile(fixture,'DEPRESION.txt'),initial);
write_matrix(fullfile(fixture,'DEPS.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'RAINFALL.txt'),zeros(3,1));
write_matrix(fullfile(fixture,'chnwidth.txt'),[1 1 0.04 1]);
write_matrix(fullfile(fixture,'link.txt'),[0 0 0 0; 0 0 0 0; 2 2 2 -1; 1 2 3 -1]);
text = setkey(base,'MODCHN','1');
text = setkey(text,'N_CHN','1');
text = setkey(text,'MAX_CHN','4');
text = setkey(text,'DEPS','DEPS.txt');
result = run_case(text,fixture);
assert(abs(sum(result.ChannelDepth(2,1:3))*2-0.4) < 1e-10);
assert(all(result.ChannelDepth(2,1:3) >= 0));

% SWMM interval construction uses seconds, constant interval inflows,
% separate hotstart files, and preserves unrelated sections.
raw = strings(24,5);
raw(1,1) = "[OPTIONS]";
raw(2,1:2) = ["FLOW_UNITS" "CMS"];
keys = ["START_DATE" "START_TIME" "END_DATE" "END_TIME" ...
    "REPORT_START_DATE" "REPORT_START_TIME" "REPORT_STEP"];
raw(3:9,1) = keys.';
raw(10,1) = "[JUNCTIONS]"; raw(11,1:3) = ["J1" "10" "2"];
raw(12,1) = "[INFLOWS]"; raw(13,1:3) = ["J1" "FLOW" "OLD"];
raw(14,1) = "[FILES]"; raw(15,1:3) = ["USE" "HOTSTART" "oldfile"];
raw(16,1) = "[TIMESERIES]"; raw(17,1:3) = ["OLD" "0:00" "1"];
swmmFile = fullfile(fixture,'swmm.inp');
write_swmm_interval(swmmFile,raw,1,0.123,datetime(2016,1,1,0,0,1),2,2);
written = fileread(swmmFile);
assert(contains(written,'00:00:03'));
assert(contains(written,'USE HOTSTART PCUW1D_HS_1'));
assert(contains(written,'SAVE HOTSTART PCUW1D_HS_0'));
assert(~contains(written,'oldfile'));
assert(contains(written,'J1 FLOW OFM_NODE_1 FLOW 1 1 0'));

% A short run on the supplied disconnected-water-body benchmark.
benchmark = fileread(fullfile(root,'Test_01.inf'));
benchmark = setkey(benchmark,'END_TIME','2016 7 5 0 0 10');
benchmark = setkey(benchmark,'DELTA_T','1');
benchmark = setkey(benchmark,'RES_GEN',['regression_' char(java.util.UUID.randomUUID)]);
result = run_case(benchmark,fixture);
assert(all(isfinite(result.Exp_flo),'all') && all(result.Exp_flo(result.Exp_flo ~= -9999) >= 0));

fprintf('All OFM regression checks passed.\n');
end

function text = setkey(text,key,value)
text = regexprep(text,['(?m)^' key '[ \t]+[^\r\n]*'],[key ' ' value]);
end

function write_matrix(file,value)
writematrix(value,file,'Delimiter','tab');
end

function result = run_case(text,fixture)
config = fullfile(fixture,'case.inf');
fid = fopen(config,'w'); fprintf(fid,'%s',text); fclose(fid);
result = ofm_urban(config);
cleanup = onCleanup(@() rmdir(fileparts(result.OutputDirectory),'s'));
end
