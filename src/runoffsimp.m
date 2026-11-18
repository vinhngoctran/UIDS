function [runoff,baseflow,slowflow,fastflow,rstorage,lstorage,perc] = ...
    runoffsimp(pa_runoff,runof_reg,infilrate,rstorage,lstorage,surfrunoff,deltaT)
%RUNOFFSIMP Reservoir runoff. Storages [m], rates [m/s], recession times [d].
% The optional step defaults to one second for older callers.
if nargin < 7, deltaT = 1; end
p = pa_runoff(runof_reg,:);
if any(~isfinite(p)) || any(p(1:4) <= 0) || p(5) < 0 || p(6) < 0
    error('OFM:RunoffInput','Invalid runoff reservoir parameters.');
end
lstorage = max(0,lstorage)+max(0,infilrate)*deltaT;
fastDepth = min(lstorage,max(0,lstorage-p(6))*deltaT/(p(1)*86400));
lstorage = lstorage-fastDepth;
slowDepth = min(lstorage,lstorage^(1+p(5))*deltaT/(p(2)*86400));
lstorage = lstorage-slowDepth;
percDepth = min(lstorage,lstorage*deltaT/(p(4)*86400));
lstorage = lstorage-percDepth;
rstorage = max(0,rstorage)+percDepth;
baseDepth = min(rstorage,rstorage*deltaT/(p(3)*86400));
rstorage = rstorage-baseDepth;
baseflow = baseDepth/deltaT;
slowflow = slowDepth/deltaT;
fastflow = fastDepth/deltaT;
perc = percDepth/deltaT;
runoff = surfrunoff+baseflow+slowflow+fastflow;
end
