function build_ofm
%BUILD_OFM Compile the threaded adaptive solver using the local C++ compiler.
% No OpenMP installation or Parallel Computing Toolbox is required.
source=fileparts(mfilename('fullpath')); root=fileparts(source);
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
old=pwd; restore=onCleanup(@() cd(old)); cd(folder);
% Use the C Matrix API entry point even with C++ compiler configurations.
flags={'-R2018a','-O','LINKEXPORTCPP='};
if isunix && ~ismac
    flags=[flags,{'CXXFLAGS=$CXXFLAGS -pthread','LDFLAGS=$LDFLAGS -pthread'}];
end
mex(flags{:},fullfile(root,'native','ofm_route.cpp'), ...
    '-outdir',source,'-output','ofm_route');
rehash;
fprintf('Built %s\n',fullfile(source,['ofm_route.' mexext]));
end
