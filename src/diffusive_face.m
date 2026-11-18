function [rate,conductance] = diffusive_face(h1,h2,z1,z2,n1,n2,d1,d2,dx,varargin)
%DIFFUSIVE_FACE Signed cell-depth transfer across a shared face [m/s].
% Ordinary faces use Manning's law. A roof/ground face uses the original
% UIDS roof outlet law Q = 0.6*0.5*sqrt(2*g)*h_roof^(3/2) [m^3/s].
% Optional logical roof masks distinguish buildings from terrain elevation.
% Without masks, DEM=999 retains the original building-marker convention.
if isempty(varargin)
    roof1 = z1 == 999;
    roof2 = z2 == 999;
else
    if numel(varargin) < 2
        error('OFM:RoofMask','Pass both face roof masks.');
    end
    roof1 = logical(varargin{1});
    roof2 = logical(varargin{2});
end
mixed = xor(roof1,roof2);
head = (z1+h1)-(z2+h2);
forward = head >= 0;
h = h1; n = n1; dep = d1;
h(~forward) = h2(~forward);
n(~forward) = n2(~forward);
dep(~forward) = d2(~forward);
faceDepth = max(0,min(h-dep,max(z1+h1,z2+h2)-max(z1,z2)));
rate = zeros(size(head));
conductance = rate;
wet = faceDepth > 0 & abs(head) > 1e-12 & ~mixed;
if any(~isfinite(n(wet)) | n(wet) <= 0)
    error('OFM:Roughness','Manning coefficients must be finite and positive.');
end
% Divide the volume flux by cell area to obtain a depth-transfer rate.
rate(wet) = sign(head(wet)).*faceDepth(wet).^(5/3).* ...
    sqrt(abs(head(wet))/dx)./(n(wet)*dx);
conductance(wet) = abs(rate(wet))./abs(head(wet)) + ...
    (5/3)*abs(rate(wet))./faceDepth(wet);

% Preserve the original sequential scheme's steep-slope outlet law.
% The adaptive diffusion solver uses Manning conveyance for ordinary faces.
if numel(varargin) >= 3 && varargin{3}
    steep = wet & abs(head)/dx >= 0.1;
    rate(steep) = sign(head(steep)).*0.6.*faceDepth(steep).^(3/2)/dx;
    conductance(steep) = abs(rate(steep))./abs(head(steep)) + ...
        1.5*abs(rate(steep))./faceDepth(steep);
end

% Roof drainage is driven by water above roof depression storage, not by
% the artificial 999-metre drop. Ground water does not enter through this
% simplified outlet. Each adjacent roof/ground face represents an outlet,
% as in the original UIDS implementation; no downpipe network is inferred.
roofDepth = h1-d1;
roofDepth(~roof1) = h2(~roof1)-d2(~roof1);
roofDepth = max(0,roofDepth);
outlet = mixed & roofDepth > 0;
coefficient = 0.6*0.5*sqrt(2*9.81)/dx^2;
direction = 2*double(roof1)-1;
rate(outlet) = direction(outlet).*coefficient.*roofDepth(outlet).^(3/2);
conductance(outlet) = 1.5*coefficient*sqrt(roofDepth(outlet));
end
