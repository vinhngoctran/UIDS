function build_swmm
%BUILD_SWMM Compile the pinned EPA engine and in-memory MATLAB interface.
% Run once on each platform with a supported C compiler (mex -setup C).
source = fileparts(mfilename('fullpath'));
root = fileparts(source);
engine = fullfile(root,'third_party','swmm');
files = dir(fullfile(engine,'*.c'));
files(strcmp({files.name},'swmm5.c')) = []; % included by the adapter
units = arrayfun(@(f) fullfile(f.folder,f.name),files,'UniformOutput',false);
buildFolder = tempname;
mkdir(buildFolder);
cleanup = onCleanup(@() rmdir(buildFolder,'s'));
old = pwd;
restore = onCleanup(@() cd(old));
cd(buildFolder);
mex('-R2018a','-O',['-I' engine],['-I' fullfile(engine,'include')], ...
    fullfile(root,'native','swmm_live.c'),units{:}, ...
    '-outdir',source,'-output','swmm_live');
rehash;
fprintf('Built %s\n',fullfile(source,['swmm_live.' mexext]));
end
