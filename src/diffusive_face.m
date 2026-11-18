function [rate, conductance] = diffusive_face(h1,h2,z1,z2,n1,n2,d1,d2,dx)
%DIFFUSIVE_FACE Signed cell-depth transfer rate across a shared face [m/s].
% Positive flow goes from cell 1 to cell 2. Manning's law uses water
% surface slope. The higher bed blocks flow until water overtops it.
head = (z1 + h1) - (z2 + h2);
forward = head >= 0;
h = h1; n = n1; dep = d1;
h(~forward) = h2(~forward);
n(~forward) = n2(~forward);
dep(~forward) = d2(~forward);
faceDepth = max(0, min(h-dep, max(z1+h1,z2+h2)-max(z1,z2)));
rate = zeros(size(head));
conductance = rate;
wet = faceDepth > 0 & abs(head) > 1e-12;
if any(~isfinite(n(wet)) | n(wet) <= 0)
    error('OFM:Roughness','Manning coefficients must be finite and positive.');
end
% q = dx*h^(5/3)*sqrt(abs(deta/dx))/n [m^3/s]. Divide by dx^2.
rate(wet) = sign(head(wet)).*faceDepth(wet).^(5/3).* ...
    sqrt(abs(head(wet))/dx)./(n(wet)*dx);
% A diffusion bound plus the derivative of the depth-dependent conveyance.
conductance(wet) = abs(rate(wet))./abs(head(wet)) + ...
    (5/3)*abs(rate(wet))./faceDepth(wet);
end
