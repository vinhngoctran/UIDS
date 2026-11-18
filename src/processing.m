function data = processing(dataraw,deltaT,dt,t_dur,type)
%PROCESSING Resample native records. Types: 0 mean, 1 total, 2 linear.
validateattributes(deltaT,{'numeric'},{'scalar','positive','finite'});
steps = round(t_dur/deltaT);
if abs(steps*deltaT-t_dur) > 1e-8*max(1,t_dur)
    error('OFM:TimeStep','Duration must be a whole number of DELTA_T steps.');
end
gridded = ndims(dataraw) == 3;
if gridded
    data = zeros(size(dataraw,1),size(dataraw,2),steps);
else
    data = zeros(steps,size(dataraw,2));
end
for k = 1:steps
    value = sample_forcing(dataraw,k,deltaT,dt,type,gridded);
    if gridded, data(:,:,k) = value; else, data(k,:) = value; end
end
end