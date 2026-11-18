% OFM-Urban: run one selected configuration and write each step to disk.
% Start this script from the repository root.
clearvars; close all; clc;
addpath(fullfile(pwd,'src'));
Filename = 'Test_01.inf';
Results = ofm_urban(Filename,struct('Plot',true));
fprintf('Step snapshots and summary: %s\n',Results.OutputDirectory);

% Extract point series without loading all gridded results into memory.
% For Test_01, sample every 60 seconds (NT_REP in the original example).
metadata = load(fullfile(Results.OutputDirectory,'metadata.mat'),'deltaT');
stride = max(1,round(60/metadata.deltaT));
[timeSeconds,waterDepth] = read_ofm_series(Results.OutputDirectory, ...
    [7 43; 7 63],stride:stride:Results.StepsWritten);
figure;
plot(timeSeconds/3600,waterDepth,'LineWidth',2);
xlabel('Time [hours]'); ylabel('Water depth [m]');
legend('Cell (7,43)','Cell (7,63)');

% Other configurations: Test_02.inf, Test_03.inf, Test_05.inf,
% Test_08A.inf, GN2010.inf, GN2011.inf, GN2013.inf.
% Test_03 includes momentum effects outside the diffusive-wave approximation.
