function rate = overland(TSP2D,deltaT,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,of_rate,varargin)
%OVERLAND Legacy sequential routing using the predicted cell depths.
% TSTP2D=1 is a sequential, availability-limited update, not a globally
% implicit diffusion solve. Use TSTP2D=2 for adaptive finite-volume routing.
h1 = max(0,oflow(x,y)+of_rate(x,y)*deltaT);
h2 = max(0,oflow(xx,yy)+of_rate(xx,yy)*deltaT);
rate = diffusive_face(h1,h2,dem(x,y),dem(xx,yy),man_area(x,y), ...
    man_area(xx,yy),s_dep(x,y),s_dep(xx,yy),s_grid);
if rate >= 0
    available = max(0,h1-s_dep(x,y));
else
    available = max(0,h2-s_dep(xx,yy));
end
rate = sign(rate)*min(abs(rate),available/deltaT);
end
