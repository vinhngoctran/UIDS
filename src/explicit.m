function [of_rate,uvfluxs] = explicit(TSP2D,shp,n_row,n_col,deltaT,s_dep,s_grid,man_area,dem,oflow,of_rate,varargin)
%EXPLICIT Conservative finite-volume diffusive-wave routing.
% Rates, including the input source/sink rate, are depths per second.
% Sources are applied once, then routing uses adaptive explicit substeps.
% uvfluxs(:,:,1:2) are mean east/south depth transfer rates [m/s].
validateattributes(deltaT,{'numeric'},{'scalar','positive','finite'});
validateattributes(s_grid,{'numeric'},{'scalar','positive','finite'});
active = shp ~= -9999;
roofs = active & (shp == 999 | dem == 999);
fixed = active & (shp == 4 | shp == 5);
fixedDepth = oflow;
fixedDepth(shp == 5) = 0;
if ~isempty(varargin)
    fixedDepth = varargin{1};
end
% Inactive cells must not enter face calculations with no-data roughness.
dem(~active) = 0; man_area(~active) = 1;
s_dep(~active) = 0; oflow(~active) = 0;
oldDepth = oflow;
if any(~isfinite(oflow(active))) || any(~isfinite(of_rate(active)))
    error('OFM:NonfiniteState','Surface depth and source rates must be finite.');
end
oflow(active) = max(0,oflow(active) + deltaT*of_rate(active));
oflow(fixed) = fixedDepth(fixed);
remaining = deltaT;
uvfluxs = zeros(n_row,n_col,2);
east = active(:,1:end-1) & active(:,2:end);
south = active(1:end-1,:) & active(2:end,:);
substeps = 0;
while remaining > 0
    substeps = substeps + 1;
    if substeps > 100000
        error('OFM:Substeps','Routing exceeded 100000 substeps. Reduce DELTA_T or inspect the input.');
    end
    [qe,ce] = diffusive_face(oflow(:,1:end-1),oflow(:,2:end), ...
        dem(:,1:end-1),dem(:,2:end),man_area(:,1:end-1),man_area(:,2:end), ...
        s_dep(:,1:end-1),s_dep(:,2:end),s_grid,roofs(:,1:end-1),roofs(:,2:end));
    [qs,cs] = diffusive_face(oflow(1:end-1,:),oflow(2:end,:), ...
        dem(1:end-1,:),dem(2:end,:),man_area(1:end-1,:),man_area(2:end,:), ...
        s_dep(1:end-1,:),s_dep(2:end,:),s_grid,roofs(1:end-1,:),roofs(2:end,:));
    qe(~east) = 0; ce(~east) = 0;
    qs(~south) = 0; cs(~south) = 0;
    totalC = zeros(n_row,n_col);
    outgoing = totalC;
    totalC(:,1:end-1) = totalC(:,1:end-1) + ce;
    totalC(:,2:end) = totalC(:,2:end) + ce;
    totalC(1:end-1,:) = totalC(1:end-1,:) + cs;
    totalC(2:end,:) = totalC(2:end,:) + cs;
    outgoing(:,1:end-1) = outgoing(:,1:end-1) + max(qe,0);
    outgoing(:,2:end) = outgoing(:,2:end) + max(-qe,0);
    outgoing(1:end-1,:) = outgoing(1:end-1,:) + max(qs,0);
    outgoing(2:end,:) = outgoing(2:end,:) + max(-qs,0);
    dt = remaining;
    maxC = max(totalC(:));
    if maxC > 0, dt = min(dt,0.4/maxC); end
    draining = active & ~fixed & outgoing > 0;
    available = max(0,oflow-s_dep);
    if any(draining(:))
        dt = min(dt,0.9*min(available(draining)./outgoing(draining)));
    end
    if ~isfinite(dt) || dt <= 0 || remaining-dt == remaining
        error('OFM:SubstepUnderflow','Routing cannot advance. Inspect depths, roughness and grid size.');
    end
    change = zeros(n_row,n_col);
    change(:,1:end-1) = change(:,1:end-1) - qe;
    change(:,2:end) = change(:,2:end) + qe;
    change(1:end-1,:) = change(1:end-1,:) - qs;
    change(2:end,:) = change(2:end,:) + qs;
    oflow(active) = oflow(active) + dt*change(active);
    oflow(fixed) = fixedDepth(fixed);
    uvfluxs(:,1:end-1,1) = uvfluxs(:,1:end-1,1) + dt*qe/deltaT;
    uvfluxs(1:end-1,:,2) = uvfluxs(1:end-1,:,2) + dt*qs/deltaT;
    remaining = max(0,remaining-dt);
end
of_rate = zeros(n_row,n_col);
of_rate(active) = (oflow(active)-oldDepth(active))/deltaT;
end
