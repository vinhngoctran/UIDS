%============================================================================
%   explicit.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate overland flow by 2D diffusive and continuity equations
%       with explicit scheme for adaptive timestepping
%============================================================================

function [of_rate,uvfluxs] = explicit(TSP2D,shp,n_row,n_col,deltaT,s_dep,s_grid,man_area,dem,oflow,of_rate)

% initialize
old_oflow = oflow;
next_oflow = oflow;
dt_next = 0;
max_vec = s_grid/deltaT;
t = 0;
while (t < deltaT)
    new_of_rate = of_rate;
    dt = deltaT - t;
    reset_dt(1:n_row,1:n_col) = dt;
    
    % Compute the flow rate (u,v) using 2D diffusive and continuity equations
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) ~= -9999
                for z = -1:0
                    xx = x + z + 1;
                    yy = y - z;
                    
                    %Compute new delta T using a coefficient of 0.7 (this coefficient can vary between 0 and 1)
                     % Older
                    if max_vec * dt > s_grid
                        dt = 0.5 * s_grid / max_vec;  % apha = 0.5 == stability coefficients
                    end
                     % New
%                     if oflow(x,y) > 0
%                         dt = 0.5* s_grid / sqrt(9.81 * oflow(x,y)); % New
%                     end
                    if dt_next ~= 0
                        dt = dt_next;
                    end
                    if dt < 0.001
                        dt = 0.001;
                    end
                    
%                     if max_vec * dt > s_grid
%                         pdt = 0.5* s_grid / sqrt(9.81 * oflow(x,y));
%                     elseif dt_next~= 0
%                         pdt = dt_next;
%                     else
%                         pdt = dt;
%                     end
%                     
%                     pdt = max([0.01,pdt]);
%                     
%                     if pdt < dt
%                         dt = pdt;
%                     end
%                     
                    if xx <= n_row && yy <= n_col && shp(xx,yy) ~= -9999
                        uvfluxs(x,y,z+2) = overland(TSP2D,dt,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,new_of_rate);
                        new_of_rate(x,y) = new_of_rate(x,y) - uvfluxs(x,y,z+2);
                        new_of_rate(xx,yy) = new_of_rate(xx,yy) + uvfluxs(x,y,z+2);
                        if max(abs(new_of_rate(x,y)),abs(new_of_rate(xx,yy))) > max_vec
                            max_vec = max(abs(new_of_rate(x,y)),abs(new_of_rate(xx,yy)));
                        end
                    end
                end
            end
        end
    end
    
    % Update water depth at time minimum dt
    for x = 1:n_row
        for y = 1:n_col
            if shp(x,y) ~= -9999
                oflow(x,y) = oflow(x,y) + new_of_rate(x,y) * dt;        % Water depth
                % Re-compute dt to avoid negative water depth
                if oflow(x,y)<0
                    if next_oflow(x,y) >  s_dep(x,y)
                        reset_dt(x,y) = 0.5* -(next_oflow(x,y) - s_dep(x,y))/ (new_of_rate(x,y));
                    else
                        reset_dt(x,y) = 0.5 * (s_dep(x,y) - next_oflow(x,y))/ (of_rate(x,y));
                    end
                end
            end
        end
    end
%     dt
%     t
%     min(oflow(:))
    % Return to the previous time step
    if min(oflow(:)) < 0 && dt > 0.001
        [aa bb] = find(oflow == min(oflow(:)));
        oflow = next_oflow;
        dt_next = min(reset_dt(:));
    elseif min(oflow(:)) > -0.0001
        next_oflow = oflow;
        dt_next = 0;
        % update timestep
        t = t + dt;
    end
% dt_next
end

% Normalize flow rate (of_rate) according to deltaT
for x = 1:n_row
    for y = 1:n_col
        if shp(x,y) ~= -9999
            of_rate(x,y) = (oflow(x,y) - old_oflow(x,y)) / deltaT;        % Water depth
        end
    end
end
end