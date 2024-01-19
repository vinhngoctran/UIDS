%============================================================================
%   rainfest.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimates rainfall weight for grid x,y
%   Program estimate the grid rainfall from gauge measurement OR grid
%   rainfall
%============================================================================

function rain_int = rainfest(deltaT,s_grid,n_row,x,y,x_cor,y_cor,id_rain,n_rain,x_rain,y_rain,RAIN_at_t,varargin)

x_node_c = x_cor;
y_node_c = y_cor + n_row * s_grid;

if id_rain == 0           % For distributed gridded rainfall
    rain_int = RAIN_at_t(x,y);
    
else                        % For local rainfall from gauge
    
    %If there is only one rain gage, then the rainfall intensity at cell j,k 'rint(j,k)' is set equal to the rainfall to
    %the rainfall intensity for that single gage
    if n_rain == 1
        rain_int = RAIN_at_t(1);
    else
        % If there are more than one raingage (n_rain>1) then the distance from
        % cell(x,y) to raingage is computed below as 'dista'
        closest_gauge = 1;        
        xc = x_node_c + y * s_grid - s_grid/2;
        yc = y_node_c - x * s_grid + s_grid/2;
        
        % Compute distance between grid cell and rain gauge
        dist_gau = sqrt((xc - x_rain).^2 + (yc-y_rain).^2);
        
        % Find gage closet to the cell (Thiessen polygon)
        idx = find(dist_gau == min(dist_gau));
        
        rain_int = RAIN_at_t(idx);
    end
end

% Unit change from mm/deltaT to M/S
rain_int = (rain_int/1000) * 1/deltaT;

end