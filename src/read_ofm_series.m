function [timeSeconds,values] = read_ofm_series(outputDirectory,cells,steps)
%READ_OFM_SERIES Extract selected [row column] cells one snapshot at a time.
% Optional steps selects computational step numbers, for example 1:100:10000.
validateattributes(cells,{'numeric'},{'2d','ncols',2,'integer','positive'});
if nargin < 3
    files = dir(fullfile(outputDirectory,'step_*.mat'));
    [~,order] = sort({files.name});
    files = files(order);
    names = {files.name};
else
    validateattributes(steps,{'numeric'},{'vector','integer','positive'});
    names = arrayfun(@(k) sprintf('step_%09d.mat',k),steps,'UniformOutput',false);
end
timeSeconds = zeros(numel(names),1);
values = zeros(numel(names),size(cells,1));
for k = 1:numel(names)
    snapshot = load(fullfile(outputDirectory,names{k}),'Results');
    field = snapshot.Results.Exp_flo;
    idx = sub2ind(size(field),cells(:,1),cells(:,2));
    values(k,:) = field(idx).';
    timeSeconds(k) = snapshot.Results.TimeSeconds;
end
end
