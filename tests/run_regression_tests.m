function run_regression_tests
%RUN_REGRESSION_TESTS Numerical invariants and short model/output regressions.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
old = pwd; cleanup = onCleanup(@() cd(old)); cd(root);

% Closed-domain volume, symmetry, sources, depression storage and no-data.
shp = zeros(5); z = zeros(5); n = 0.04*ones(5); dep = zeros(5);
h = zeros(5); h(3,3) = 0.3;
[r,f] = explicit(2,shp,5,5,2,dep,2,n,z,h,zeros(5));
next = h+2*r;
assert(all(next(:) >= 0));
assert(abs(sum(next(:))-sum(h(:))) < 1e-12);
assert(max(abs(next-flipud(next)),[],'all') < 1e-12);
assert(max(abs(next-fliplr(next)),[],'all') < 1e-12);
assert(isequal(size(f),[5 5 2]));
[r,~] = explicit(2,shp,5,5,2,dep,2,n,z,zeros(5),ones(5)*1e-4);
assert(max(abs(r(:)-1e-4)) < 1e-12);
dep(:) = 0.05; h(:) = 0; h(3,3) = 0.3;
[r,~] = explicit(2,shp,5,5,20,dep,2,n,z,h,zeros(5));
next = h+20*r;
assert(next(3,3) >= 0.05-1e-12);
assert(abs(sum(next(:))-0.3) < 1e-12);
z = [0 2]; h = [0.5 0]; dep = [0 0]; n = [0.03 0.03];
assert(diffusive_face(h(1),h(2),z(1),z(2),n(1),n(2),dep(1),dep(2),1) == 0);
z = [1200 1200.2]; h = [0.4 0.2];
assert(abs(diffusive_face(h(1),h(2),z(1),z(2),n(1),n(2),dep(1),dep(2),1)) < 1e-10);
h = [0.4 0.1];
q1 = diffusive_face(h(1),h(2),z(1),z(2),n(1),n(2),dep(1),dep(2),1);
q2 = diffusive_face(h(2),h(1),z(2),z(1),n(2),n(1),dep(2),dep(1),1);
assert(abs(q1+q2) < 1e-12 && q1 > 0);
shp = [0 -9999]; n(2) = -9999;
[r,~] = explicit(2,shp,1,2,1,dep,1,n,z,h,[0 0]);
assert(all(r == 0));

% Original roof outlet law, conserved transfers and no ground-to-roof flow.
roofDepth = 0.12; depression = 0.02; dx = 2;
expectedQ = 0.6*0.5*sqrt(2*9.81)*(roofDepth-depression)^1.5;
r = diffusive_face(roofDepth,0,999,10,0.03,0.03,depression,0,dx);
assert(abs(r*dx^2-expectedQ) < 1e-12);
rReverse = diffusive_face(0,roofDepth,10,999,0.03,0.03,0,depression,dx);
assert(abs(r+rReverse) < 1e-12);
assert(diffusive_face(0,5,999,10,0.03,0.03,0,0,dx) == 0);
assert(diffusive_face(0.1,0.1,999,999,0.03,0.03,0,0,dx) == 0);
% Explicit masks support real roof elevations and disambiguate 999-m terrain.
rMarked = diffusive_face(roofDepth,0,15,10,0.03,0.03,depression,0,dx,true,false);
assert(abs(rMarked-r) < 1e-12);
rTerrain = diffusive_face(roofDepth,0,999,10,0.03,0.03,depression,0,dx,false,false);
assert(abs(rTerrain-r) > 1e-6);
for marked = [false true]
    shape = zeros(3); bed = zeros(3); depth = zeros(3); storage = zeros(3);
    if marked, shape(2,2) = 999; bed(2,2) = 15; else, bed(2,2) = 999; end
    depth(2,2) = 0.12; storage(2,2) = 0.02;
    [rate,~] = explicit(2,shape,3,3,60,storage,2,0.03*ones(3),bed,depth,zeros(3));
    next = depth+60*rate;
    assert(abs(sum(next,'all')-0.12) < 1e-12);
    assert(next(2,2) >= 0.02-1e-12 && next(2,2) < depth(2,2));
    assert(all(next(:) >= -1e-12));
    assert(max(abs(next-flipud(next)),[],'all') < 1e-12);
    assert(max(abs(next-fliplr(next)),[],'all') < 1e-12);
end
% Both sequential and adaptive paths use the roof outlet law.
rSequential = overland(1,0.1,[0.02 0],2,[0.03 0.03],[999 10],[0.12 0],1,1,1,2,[0 0]);
assert(abs(rSequential-r) < 1e-12);
rSteep = overland(1,0.1,[0 0],2,[0.03 0.03],[1 0],[0.1 0],1,1,1,2,[0 0]);
assert(abs(rSteep-0.6*0.1^1.5/2) < 1e-12);

% Distinct public preprocessing utilities remain available.
[dayLAI,monthLAI] = laivariant(4,25,datenum(2016,[1 6 12],1));
assert(isequal(dayLAI,[1;4;1]) && numel(monthLAI) == 12);
assert(max(abs(processing([6;12;3],4,6,12,1)-[4;6;8])) < 1e-12);

% Optional-model input helpers remain present and callable.
[latitude,longitude] = coordinate(500000,0,'52 N');
assert(abs(latitude) < 1e-8 && abs(longitude-129) < 1e-8);
soilInput = modswfip(fullfile(root,'Examples','Input','SOILFLUX.txt'));
assert(soilInput.nlayer == 1 && isequal(size(soilInput.par),[10 1]));

% Native forcing preserves cumulative totals with non-integer step ratios.
x = arrayfun(@(k) sample_forcing([6;12;3],k,4,6,1,false),1:3).';
assert(max(abs(x-[4;6;8])) < 1e-12);
x = arrayfun(@(k) sample_forcing([0;6;12],k,4,6,2,false),1:3).';
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
% Rainfall on a marked building reaches the adjacent ground in both schemes.
roofMask = zeros(3,4); roofMask(2,2) = 999;
roofBed = zeros(3,4); roofBed(2,2) = 999;
write_matrix(fullfile(fixture,'MASK.txt'),roofMask);
write_matrix(fullfile(fixture,'DEM.txt'),roofBed);
for scheme = [1 2]
    text = setkey(base,'TSTP2D',num2str(scheme));
    result = run_case(text,fixture);
    assert(abs(sum(result.Exp_flo,'all')-12*0.003) < 1e-11);
    assert(result.Exp_flo(2,2) < 0.003);
    assert(any(result.Exp_flo(roofMask ~= 999) > 0.003));
end
write_matrix(fullfile(fixture,'MASK.txt'),zeros(3,4));
write_matrix(fullfile(fixture,'DEM.txt'),zeros(3,4));

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

% Parse and prepare the supplied SWMM topology once at startup.
[network,prepared] = prepare_swmm_input(fullfile(root,'Examples', ...
    '_Benchmarking_Model_Data','Test8B','1D','1D.mdf'), ...
    datetime(2016,1,1),datetime(2016,1,1,0,0,10),0.1);
assert(isequal(network.NodeIDs,{'1','2'}));
assert(isequal(network.LinkIDs,{'1','2'}));
assert(isequal(network.OutfallIDs,{'3'}));
assert(~contains(prepared,'HOTSTART'));
assert(contains(prepared,'ROUTING_STEP 0.1'));

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
