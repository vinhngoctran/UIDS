function timing=run_routing_parallel_tests(gridSize)
%RUN_ROUTING_PARALLEL_TESTS MATLAB/native parity, conservation and timings.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
assert(exist('ofm_route','file')==3,'Run build_ofm before the native routing checks.');
if nargin<1, gridSize=256; end
rng(731);
for dims={[1 9],[9 1],[9 13]}
    shape=zeros(dims{1}); nrow=dims{1}(1); ncol=dims{1}(2);
    bed=.2*rand(dims{1}); rough=.03+.02*rand(dims{1});
    depression=.02*rand(dims{1}); water=.2*rand(dims{1}); source=1e-5*rand(dims{1});
    shape(2)=999; bed(2)=15;
    if numel(shape)>9, shape(5)=-9999; rough(5)=-9999; bed(5)=-9999; end
    for boundary=[false true]
        if boundary, shape(end)=4; shape(1)=5; end
        fixed=water; fixed(shape==5)=0;
        [reference,rf]=explicit(2,shape,nrow,ncol,.001,depression,2,rough,bed,water,source,fixed);
        for workers=[1 2 4]
            [native,nf]=ofm_route(shape,.001,depression,2,rough,bed,water,source,fixed,double(workers));
            assert(max(abs(native-reference),[],'all')<1e-11);
            assert(max(abs(nf-rf),[],'all')<1e-11);
            next=water+native*.001; active=shape~=-9999;
            assert(all(next(active)>=-1e-12));
            if ~boundary
                assert(abs(sum(next(active))-sum(water(active)+source(active)*.001))<1e-10);
            else
                assert(next(1)==0 && abs(next(end)-fixed(end))<1e-12);
            end
        end
    end
end
shape=zeros(3); shape(2,2)=999; bed=zeros(3); bed(2,2)=15;
water=zeros(3); water(2,2)=.12; depression=zeros(3); depression(2,2)=.02;
[reference,rf]=explicit(2,shape,3,3,60,depression,2,.04+zeros(3),bed,water,zeros(3),water);
for workers=[1 4]
    [native,nf]=ofm_route(shape,60,depression,2,.04+zeros(3),bed,water,zeros(3),water,double(workers));
    assert(max(abs(native-reference),[],'all')<1e-11 && max(abs(nf-rf),[],'all')<1e-11);
    next=water+60*native;
    assert(abs(sum(next,'all')-.12)<1e-11 && next(2,2)>=.02-1e-12 && next(2,2)<.12);
end
% Flat water at rest and dry obstructions must produce no transfers.
shape=zeros(3); water=.1+shape;
r=ofm_route(shape,1,shape,2,.04+shape,shape,water,shape,water,4);
assert(all(r==0,'all'));
r=ofm_route([0 0],1,[0 0],1,[.04 .04],[0 2],[.5 0],[0 0],[.5 0],4);
assert(all(r==0));
% Reject bad arguments without leaving worker threads or state unusable.
try
    bad=ofm_route(shape,1,shape,2,-1+shape,shape,water,shape,water,4); %#ok<NASGU>
    error('OFM:Test','Bad roughness was accepted.');
catch err
    assert(strcmp(err.identifier,'OFM:NativeRoutingInput'));
end
r=ofm_route(shape,1,shape,2,.04+shape,shape,water,shape,water,2);
assert(all(r==0,'all'));
% A larger reproducible field times routing only, after warming all backends.
n=gridSize; shape=zeros(n); [row,col]=ndgrid(1:n);
bed=.01*(row+col); depression=.005+shape; rough=.04+shape;
water=.1+.01*sin(row/7).*cos(col/9); source=1e-6+shape;
reference=@() explicit(2,shape,n,n,1,depression,2,rough,bed,water,source,water);
one=@() ofm_route(shape,1,depression,2,rough,bed,water,source,water,1);
many=@() ofm_route(shape,1,depression,2,rough,bed,water,source,water,4);
value=reference(); value=one(); value=many(); %#ok<NASGU>
timing=zeros(1,3);
for k=1:3
    timer=tic;
    if k==1, for j=1:5, value=reference(); end
    elseif k==2, for j=1:5, value=one(); end
    else, for j=1:5, value=many(); end
    end
    timing(k)=toc(timer)/5;
end
fprintf('ROUTING TIMING %dx%d MATLAB %.5f s, native/1 %.5f s, native/4 %.5f s\n',n,n,timing);
fprintf('All MATLAB/native routing and conservation checks passed.\n');
clear ofm_route % Also exercises orderly thread shutdown/reload.
end
