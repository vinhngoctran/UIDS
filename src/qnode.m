function q = qnode(deltaT,of_rate,s_grid,deps,oflow)
%QNODE Free inlet discharge [m^3/s], capped by available surface water.
h = max(0,oflow+of_rate*deltaT-deps);
q = min(s_grid*sqrt(2*9.81)*h^1.5,h*s_grid^2/deltaT);
end
