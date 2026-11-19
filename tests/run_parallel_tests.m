function timings=run_parallel_tests(gridSize,poolType)
%RUN_PARALLEL_TESTS Compare every saved numerical field across execution modes.
% Optional gridSize enables a larger timing fixture; defaults to a short check.
if nargin<1, gridSize=8; end
if nargin<2, poolType='Processes'; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
old=pwd; restore=onCleanup(@() cd(old)); cd(root);
assert(~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox'), ...
    'Parallel checks require Parallel Computing Toolbox.');
pool=gcp('nocreate'); owned=isempty(pool);
if owned, pool=parpool(poolType,4); end
cleanupPool=onCleanup(@() close_owned_pool(pool,owned));
fixture=tempname; mkdir(fixture); cleanup=onCleanup(@() rmdir(fixture,'s'));
base=fileread(fullfile(root,'Test_01.inf'));
n=gridSize; blank=zeros(n); [row,col]=ndgrid(1:n);
mask=blank; mask(1,1)=-9999; mask(2,2)=999;
write_matrix(fixture,'MASK.txt',mask); write_matrix(fixture,'DEM.txt',blank);
write_matrix(fixture,'DEPRESION.txt',blank); write_matrix(fixture,'IMPERVIOUS.txt',.2*ones(n));
write_matrix(fixture,'TestLAI.txt',1+row/n); write_matrix(fixture,'VegHeight.txt',.5*ones(n));
write_matrix(fixture,'tcan.txt',20+col/(5*n));
for name=["albsoil.txt","albveg.txt"], write_matrix(fixture,name,.2*ones(n)); end
for name=["emisssoil.txt","emissveg.txt"], write_matrix(fixture,name,.98*ones(n)); end
write_matrix(fixture,'RAINFALL.txt',[1;2;0]);
write_text(fullfile(fixture,'HEADING.txt'),sprintf( ...
    'UTM 52 S\nncols %d\nnrows %d\nxllcorner 500000\nyllcorner 4000000\ncellsize 2\nnodata -9999\n',n,n));
write_text(fullfile(fixture,'weather.mdf'),sprintf('TA RH SR LR WD VP PA\n20 .008 500 320 2 -9999 101.3\n'));
soilText=fileread(fullfile(root,'Examples','Input','SOILFLUX.txt'));
soilText=setkey(soilText,'nlayer','2'); soilText=setkey(soilText,'thick','300 400');
soilText=setkey(soilText,'psim_init','-33 -100'); write_text(fullfile(fixture,'SOILFLUX.txt'),soilText);
values={'FOLDER',[fixture filesep];'START_TIME','2016 7 5 12 0 0';'END_TIME','2016 7 5 12 0 1'; ...
    'DELTA_T','.5';'DELTAT_RAIN','1';'N_BOUN','0';'RES_GEN',['parallel_test_' char(java.util.UUID.randomUUID)]};
for k=1:size(values,1), base=setkey(base,values{k,1},values{k,2}); end
% Include infiltration-only reservoirs, each ET option, and layered soil.
flags=[0 0 1 0; 1 0 0 0; 0 1 0 0; 1 1 0 0; 1 0 1 0; 0 1 1 0; 1 1 1 0; 0 1 1 1; 1 0 1 1; 1 1 1 1];
if n>16, flags=flags(end,:); end
timings=zeros(size(flags,1),2);
for test=1:size(flags,1)
    text=base;
    keys={'MODINT','MODETS','MODINF','MODSWF'};
    for k=1:4, text=setkey(text,keys{k},num2str(flags(test,k))); end
    file=fullfile(fixture,'model.inf'); write_text(file,text);
    for execution=1:2
        options=struct('Parallel',execution==2,'NumWorkers',min(4,pool.NumWorkers));
        timer=tic; result=ofm_urban(file,options); timings(test,execution)=toc(timer);
        if execution==1
            serial=result; cleanupSerial=onCleanup(@() remove_output_folder(fileparts(serial.OutputDirectory)));
        else
            cleanupParallel=onCleanup(@() remove_output_folder(fileparts(result.OutputDirectory)));
            assert(result.Execution.ParallelHydrology && result.Execution.Workers>=2);
            for step=1:serial.StepsWritten
                filename=sprintf('step_%09d.mat',step);
                a=load(fullfile(serial.OutputDirectory,filename));
                b=load(fullfile(result.OutputDirectory,filename));
                compare_numeric(a.Results,b.Results);
                assert(max(abs(a.Waterdepth-b.Waterdepth),[],'all')<1e-12);
                assert(max(abs(a.WaterVec-b.WaterVec),[],'all')<1e-12);
            end
            clear cleanupParallel cleanupSerial
        end
    end
    fprintf('PARALLEL TIMING %dx%d serial %.3f s, %s/%d workers %.3f s, ratio %.2fx\n', ...
        n,n,timings(test,1),result.Execution.PoolType,result.Execution.Workers,timings(test,2),timings(test,1)/timings(test,2));
end
fprintf('All serial/parallel snapshot comparisons passed.\n');
end

function compare_numeric(a,b)
names=fieldnames(a);
for j=1:numel(names)
    key=names{j};
    if isnumeric(a.(key)) || islogical(a.(key))
        assert(isequal(size(a.(key)),size(b.(key))),'Shape mismatch: %s',key);
        assert(all(abs(a.(key)-b.(key))<=1e-12*max(1,abs(a.(key))),'all'),'Numerical mismatch: %s',key);
    end
end
end
function write_matrix(folder,name,data)
writematrix(data,fullfile(folder,name),'Delimiter','tab');
end
function write_text(file,text)
f=fopen(file,'w'); assert(f>=0); cleanup=onCleanup(@() fclose(f)); fprintf(f,'%s',text);
end
function text=setkey(text,key,value)
text=regexprep(text,['(?m)^' key '[ \t]+[^\r\n]*'],[key ' ' value]);
end
function close_owned_pool(pool,owned)
if owned && isvalid(pool), delete(pool); end
end

function remove_output_folder(folder)
if exist(folder,'dir'), rmdir(folder,'s'); end
end
