function rate = overland(TSP2D,deltaT,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,of_rate,varargin)
%OVERLAND Legacy sequential routing using the predicted cell depths.
% TSTP2D=1 is a sequential, availability-limited update, not a globally
% implicit diffusion solve. Use TSTP2D=2 for adaptive finite-volume routing.
roof1 = dem(x,y) == 999;
roof2 = dem(xx,yy) == 999;
if ~isempty(varargin)
    roofs = varargin{1};
    roof1 = roofs(x,y); roof2 = roofs(xx,yy);
end
h1 = max(0,oflow(x,y)+of_rate(x,y)*deltaT);
h2 = max(0,oflow(xx,yy)+of_rate(xx,yy)*deltaT);
rate = diffusive_face(h1,h2,dem(x,y),dem(xx,yy),man_area(x,y), ...
    man_area(xx,yy),s_dep(x,y),s_dep(xx,yy),s_grid,roof1,roof2,TSP2D == 1);
if rate >= 0
    available = max(0,h1-s_dep(x,y));
else
    available = max(0,h2-s_dep(xx,yy));
end
rate = sign(rate)*min(abs(rate),available/deltaT);
end
