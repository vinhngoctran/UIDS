function rate = overland2(TSP2D,deltaT,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,varargin)
%OVERLAND2 Signed depth transfer rate [m/s], positive toward (xx,yy).
rate = diffusive_face(oflow(x,y),oflow(xx,yy),dem(x,y),dem(xx,yy), ...
    man_area(x,y),man_area(xx,yy),s_dep(x,y),s_dep(xx,yy),s_grid);
if TSP2D == 1
    if rate >= 0
        available = max(0,oflow(x,y)-s_dep(x,y));
    else
        available = max(0,oflow(xx,yy)-s_dep(xx,yy));
    end
    rate = sign(rate)*min(abs(rate),available/deltaT);
end
end
