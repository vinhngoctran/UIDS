function [q,of_rate] = qmanhole(deltaT,s_grid,deps,oflow,of_rate,node_depth_at_t,dem,node_invert)
%QMANHOLE Surface-to-sewer discharge [m^3/s]. node_invert is rim depth [m].
h = max(0,oflow+of_rate*deltaT-deps);
sewerHead = node_depth_at_t+dem-node_invert;
if h == 0 || sewerHead >= dem+h
    q = 0;
elseif sewerHead < dem
    q = 1.5*sqrt(2*9.81)*h^1.5;
else
    q = 1.5^2*sqrt(2*9.81)*sqrt(max(0,dem+h-sewerHead));
endq = min(q,h*s_grid^2/deltaT);
of_rate = of_rate-q/s_grid^2;
end
