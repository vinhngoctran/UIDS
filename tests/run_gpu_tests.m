function timing=run_gpu_tests(gridSize)
%RUN_GPU_TESTS Metal/CPU accuracy, conservation, complete runs and timing.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
old=pwd;restore=onCleanup(@() cd(old));cd(root);
if nargin<1, gridSize=256; end
assert(exist('ofm_metal','file')==3,'Run build_gpu before the GPU checks.');
info=ofm_metal('info'); assert(info.Available); fprintf('Metal device: %s\n',info.Name);
rng(1801);
for dims={[1 9],[9 1],[9 13]}
    sz=dims{1}; shape=zeros(sz); bed=1200+.2*rand(sz); h=.2*rand(sz);
    dep=.02*rand(sz); rough=.03+.02*rand(sz); source=1e-5*rand(sz);
    shape(2)=999;
    if numel(shape)>9, shape(5)=-9999; bed(5)=-9999; rough(5)=-9999; end
    for boundary=[false true]
        if boundary, shape(1)=5; shape(end)=4; end
        fixed=h; fixed(shape==5)=0;
        [reference,rf]=explicit(2,shape,sz(1),sz(2),.001,dep,2,rough,bed,h,source,fixed);
        [gpu,gf,d]=ofm_metal(shape,.001,dep,2,rough,bed,h,source,fixed);
        assert(max(abs((gpu-reference)*.001),[],'all')<1e-7);
        assert(max(abs(gf-rf),[],'all')<2e-6);
        active=shape~=-9999; next=h+.001*gpu;
        assert(all(next(active)>=-1e-12) && abs(d.BalanceResidual)<1e-10);
        if ~boundary
            assert(abs(sum(next(active))-sum(h(active)+source(active)*.001))<1e-10);
        else
            assert(abs(next(1))<1e-12 && abs(next(end)-fixed(end))<1e-12);
        end
    end
end
[single,~,d]=ofm_metal(0,.5,0,2,.04,0,.1,1e-3,.1);
assert(abs(single-1e-3)<1e-12 && abs(d.BalanceResidual)<1e-12);
% High elevations and a small surface gradient must not cancel in FP32.
shape=[0 0];bed=[1200 1200.2];h=[.4 .2];dep=[0 0];rough=[.03 .03];
[gpu,~,d]=ofm_metal(shape,.1,dep,1,rough,bed,h,[0 0],h);
assert(max(abs(gpu),[],'all')<1e-8 && abs(d.BalanceResidual)<1e-12);
% Dry obstruction and a dry roof cannot receive water through an outlet.
[gpu,~,d]=ofm_metal(shape,1,dep,1,rough,[0 2],[.5 0],[0 0],[.5 0]);
assert(all(gpu==0));
[gpu,~,d]=ofm_metal([999 0],1,dep,1,rough,[15 0],[0 .5],[0 0],[0 .5]);
assert(all(gpu==0));
% Long adaptive roof drainage: all four outlets share a bounded donor.
shape=zeros(3);shape(2,2)=999;bed=zeros(3);bed(2,2)=15;
h=zeros(3);h(2,2)=.12;dep=zeros(3);dep(2,2)=.02;rough=.04+zeros(3);
[reference,rf]=explicit(2,shape,3,3,60,dep,2,rough,bed,h,zeros(3),h);
[gpu,gf,d]=ofm_metal(shape,60,dep,2,rough,bed,h,zeros(3),h);next=h+60*gpu;
assert(max(abs((gpu-reference)*60),[],'all')<2e-6);
assert(abs(sum(next,'all')-.12)<1e-10 && next(2,2)>=.02-1e-10 && next(2,2)<.12);
assert(abs(d.BalanceResidual)<1e-10 && d.Substeps>1);
% Invalid inputs must leave the device usable for another run.
try
    bad=ofm_metal(shape,1,dep,2,-rough,bed,h,zeros(3),h); %#ok<NASGU>
    error('OFM:Test','Invalid roughness was accepted.');
catch err
    assert(strcmp(err.identifier,'OFM:MetalInput'));
end
[gpu,~,d]=ofm_metal(zeros(3),1,zeros(3),2,rough,zeros(3),.1+zeros(3),zeros(3),.1+zeros(3));
assert(max(abs(gpu),[],'all')<1e-12 && abs(d.BalanceResidual)<1e-12);

% Complete driver: GPU routing with streamed output and optional hydrology.
fixture=tempname;mkdir(fixture);clean=onCleanup(@() rmdir(fixture,'s'));
n=6; blank=zeros(n); mask=blank;mask(1,1)=-9999;mask(2,2)=999;
write_matrix(fixture,'MASK.txt',mask);write_matrix(fixture,'DEM.txt',blank);
write_matrix(fixture,'DEPRESION.txt',blank);write_matrix(fixture,'IMPERVIOUS.txt',.2+blank);
write_matrix(fixture,'TestLAI.txt',2+blank);write_matrix(fixture,'VegHeight.txt',.5+blank);
write_matrix(fixture,'tcan.txt',20+blank);
for name=["albsoil.txt","albveg.txt"],write_matrix(fixture,name,.2+blank);end
for name=["emisssoil.txt","emissveg.txt"],write_matrix(fixture,name,.98+blank);end
write_matrix(fixture,'RAINFALL.txt',[1;1]);
write_text(fullfile(fixture,'HEADING.txt'),sprintf( ...
    'UTM 52 S\nncols %d\nnrows %d\nxllcorner 500000\nyllcorner 4000000\ncellsize 2\nnodata -9999\n',n,n));
write_text(fullfile(fixture,'weather.mdf'),sprintf('TA RH SR LR WD VP PA\n20 .008 500 320 2 -9999 101.3\n'));
soil=fileread(fullfile(root,'Examples','Input','SOILFLUX.txt'));
soil=setkey(soil,'nlayer','2');soil=setkey(soil,'thick','300 400');soil=setkey(soil,'psim_init','-33');
write_text(fullfile(fixture,'SOILFLUX.txt'),soil);
base=fileread(fullfile(root,'Test_01.inf'));
settings={'FOLDER',[fixture filesep];'START_TIME','2016 7 5 12 0 0';'END_TIME','2016 7 5 12 0 2'; ...
    'DELTA_T','.5';'DELTAT_RAIN','1';'N_BOUN','0';'RES_GEN',['gpu_test_' char(java.util.UUID.randomUUID)]};
for k=1:size(settings,1),base=setkey(base,settings{k,1},settings{k,2});end
for hydrology=[false true]
    text=base;
    for key={'MODINT','MODETS','MODINF','MODSWF'},text=setkey(text,key{1},num2str(hydrology));end
    file=fullfile(fixture,'case.inf');write_text(file,text);
    cpu=ofm_urban(file,struct('Parallel',false,'Routing','matlab'));
    outputCleanup=onCleanup(@() remove_outputs(cpu.OutputDirectory));
    gpu=ofm_urban(file,struct('Device','gpu','Parallel',false));
    assert(strcmp(gpu.Execution.Device,'gpu') && strcmp(gpu.Execution.RoutingBackend,'metal'));
    for step=1:gpu.StepsWritten
        filename=sprintf('step_%09d.mat',step);
        a=load(fullfile(cpu.OutputDirectory,filename));b=load(fullfile(gpu.OutputDirectory,filename));
        assert(max(abs(a.Waterdepth-b.Waterdepth),[],'all')<2e-6);
        assert(max(abs(a.WaterVec-b.WaterVec),[],'all')<2e-6);
        assert(abs(b.Results.RoutingDiagnostics.BalanceResidual)<1e-10);
        names=fieldnames(a.Results);
        for j=1:numel(names)
            name=names{j}; av=a.Results.(name); bv=b.Results.(name);
            if isnumeric(av) || islogical(av)
                assert(isequal(size(av),size(bv)) && all(abs(av-bv)<=2e-6*max(1,abs(av)),'all'), ...
                    'CPU/GPU state mismatch: %s',name);
            end
        end
        state=b.Results; active=mask~=-9999;
        budget=state.Exp_flo+state.total_evap;
        initial=0;
        if hydrology
            budget=budget+state.canopystore+state.snowpackstore+(state.swat+state.gwat)/1000+state.seep_total;
            initial=.3*.9*700/1000;
            assert(max(abs(state.soil_balance(active)))<1e-8);
        end
        assert(abs(sum(budget(active))-nnz(active)*(initial+step*.0005))<1e-9);
    end
    clear outputCleanup
    if hydrology && ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox')
        parallelGPU=ofm_urban(file,struct('Device','gpu','Parallel',true,'NumWorkers',4));
        parallelCleanup=onCleanup(@() remove_outputs(parallelGPU.OutputDirectory));
        assert(parallelGPU.Execution.ParallelHydrology && strcmp(parallelGPU.Execution.Device,'gpu'));
        assert(max(abs(parallelGPU.Exp_flo-gpu.Exp_flo),[],'all')<1e-12);
        assert(max(abs(parallelGPU.psim-gpu.psim),[],'all')<1e-10);
        clear parallelCleanup
    end
end

text=setkey(base,'TSTP2D','1');file=fullfile(fixture,'unsupported.inf');write_text(file,text);
try
    invalid=ofm_urban(file,struct('Device','gpu')); %#ok<NASGU>
    error('OFM:Test','GPU request accepted a sequential scheme.');
catch err
    assert(strcmp(err.identifier,'OFM:RoutingScheme'));
end

% Warm routing-only timings. No model/disk/worker-pool startup included.
n=gridSize;shape=zeros(n);[row,col]=ndgrid(1:n);bed=.01*(row+col);
h=.1+.01*sin(row/7).*cos(col/9);dep=.005+shape;rough=.04+shape;source=1e-6+shape;
cpu=@() ofm_route(shape,1,dep,2,rough,bed,h,source,h,4);
gpu=@() ofm_metal(shape,1,dep,2,rough,bed,h,source,h);
value=cpu();value=gpu(); %#ok<NASGU>
timing=zeros(1,2);
for k=1:2
    timer=tic;
    for j=1:3
        if k==1,value=cpu();else,value=gpu();end %#ok<NASGU>
    end
    timing(k)=toc(timer)/3;
end
fprintf('GPU TIMING %dx%d CPU/4 %.6f s, Metal %.6f s, ratio %.2fx\n',n,n,timing,timing(1)/timing(2));
fprintf('All Metal GPU accuracy, water-balance and driver checks passed.\n');
end
function write_matrix(folder,name,data)
writematrix(data,fullfile(folder,name),'Delimiter','tab');
end
function write_text(file,text)
f=fopen(file,'w');assert(f>=0);cleanup=onCleanup(@() fclose(f));fprintf(f,'%s',text);
end
function text=setkey(text,key,value)
text=regexprep(text,['(?m)^' key '[ \t]+[^\r\n]*'],[key ' ' value]);
end
function remove_outputs(folder)
folder=fileparts(folder);if exist(folder,'dir'),rmdir(folder,'s');end
end
