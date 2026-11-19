function build_gpu
%BUILD_GPU Build the Apple silicon Metal backend for adaptive OFM routing.
% Uses the local C++ compiler and Apple frameworks; no CUDA/toolbox needed.
if ~ismac || ~strcmp(computer('arch'),'maca64')
    error('OFM:GPUPlatform','The current GPU backend requires MATLAB on Apple silicon. CPU routing remains available on other platforms.');
end
source=fileparts(mfilename('fullpath')); root=fileparts(source);
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
shader=fileread(fullfile(root,'native','ofm_metal_kernels.metal'));
header=fullfile(folder,'ofm_metal_kernels.hpp'); f=fopen(header,'w');
if f<0, error('OFM:GPUBuild','Cannot create the embedded shader header.'); end
fileCleanup=onCleanup(@() fclose(f));
fprintf(f,'static const char *ofmMetalSource=R"OFM_METAL(%s)OFM_METAL";\n',shader); clear fileCleanup
% mex selects its C++ configuration for .cpp, with Objective-C++ enabled.
cpp=fullfile(folder,'ofm_metal.cpp'); copyfile(fullfile(root,'native','ofm_metal.mm'),cpp);
old=pwd; restore=onCleanup(@() cd(old)); cd(folder);
mex('-R2018a','-O','LINKEXPORTCPP=', ...
    'CXXFLAGS=$CXXFLAGS -x objective-c++ -fobjc-arc', ...
    'LDFLAGS=$LDFLAGS -framework Foundation -framework Metal', ...
    ['-I' folder],cpp,'-outdir',source,'-output','ofm_metal');
rehash;
fprintf('Built %s\n',fullfile(source,['ofm_metal.' mexext]));
end
