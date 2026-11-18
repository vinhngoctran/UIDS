function value = sample_forcing(raw,step,deltaT,nativeDt,kind,gridded)
%SAMPLE_FORCING Read one model interval without expanding the time series.
% kind: 0 interval mean, 1 interval accumulation, 2 end-of-step interpolation.
if nargin < 6, gridded = false; end
validateattributes(nativeDt,{'numeric'},{'scalar','positive','finite'});
a = (step-1)*deltaT/nativeDt;
b = step*deltaT/nativeDt;
if gridded
    count = size(raw,3);
else
    count = size(raw,1);
end
if kind == 2
    left = floor(b+1e-10);
    fraction = b-left;
    if abs(fraction) < 1e-10
        value = take(raw,left+1,gridded,count);
    else
        value = (1-fraction)*take(raw,left+1,gridded,count) + ...
            fraction*take(raw,left+2,gridded,count);
    end
else
    value = zeros(size(take(raw,1,gridded,count)));
    for k = floor(a+1e-10)+1:ceil(b-1e-10)
        weight = max(0,min(b,k)-max(a,k-1));
        value = value + weight*take(raw,k,gridded,count);
    end
    if kind == 0
        value = value/(b-a);
    elseif kind ~= 1
        error('OFM:ForcingType','Unknown forcing type.');
    end
end
end

function value = take(raw,k,gridded,count)
if k < 1 || k > count
    error('OFM:ForcingLength','Forcing record %d is required, but only %d records exist.',k,count);
end
if gridded, value = raw(:,:,k); else, value = raw(k,:); end
end
