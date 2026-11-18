function q = qoutlet(s_grid,deltaT,oflow,x_out,y_out,iout,of_rate,D_out,W_out,pflow)
%QOUTLET Manning outlet flow [m^3/s], capped for every positive timestep.
x = x_out(iout); y = y_out(iout);
h = max(0,oflow(x,y)+of_rate(x,y)*deltaT-D_out(iout));
q = min(W_out(iout)*pflow*h^(5/3),h*s_grid^2/deltaT);
end
