%============================================================================
%   overland.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate overland flow by 2D diffusive and continuity equations
%   with implicit scheme
%
%s_grid: THE GRID SIZE
%man_area:
%Manning: MANNING COEFFICIENT
%STEP
%S0: BED SLOPE
%slop_fic: FRICTION SLOPE
%slop_wl: WATER LEVEL SLOPE
%oflow: WATER LEVEL
%DQQ: UNIT OVERLAND FLOW RATE
%of_rate: FLOWRATE IN EACH GRID
%============================================================================

function velocity = overland2(TSP2D,deltaT,s_dep,s_grid,man_area,dem,oflow,x,y,xx,yy,varargin)

s0 = (dem(x,y) - dem(xx,yy))/s_grid;
slop_wl = ((oflow(xx,yy)) - (oflow(x,y)))/s_grid;
slop_fic = s0 - slop_wl;
pos_oflow = oflow(x,y);
rman = man_area(x,y);
deps = s_dep(x,y);
if slop_fic < 0
    pos_oflow = oflow(xx,yy);
    rman = man_area(xx,yy);
    deps = s_dep(xx,yy);
end

%% Compute OFM in ROOF
if dem(x,y) == 999 && dem(xx,yy) == 999

    if pos_oflow <= deps || slop_fic == 0
        dqq = 0;

    else
        alfa = ((abs(slop_fic))^0.5) / rman;
        h_ofm = pos_oflow - deps;
        
        if slop_fic >= 0
            dqq = s_grid * alfa * (h_ofm^(5/3));
        else
            dqq = -s_grid * alfa * (h_ofm.^(5/3));
        end

    end
    
%% Compute OFM in LAND
elseif dem(x,y) < 999 && dem(xx,yy) < 999  
    
    if pos_oflow <= deps || slop_fic == 0
        dqq = 0;        
    else
        alfa = ((abs(slop_fic))^0.5) / rman;
        h_ofm = pos_oflow - deps;
        
        if slop_fic >= 0
            dqq = s_grid * alfa * (h_ofm^(5/3));
        else
            dqq = -s_grid * alfa * (h_ofm.^(5/3));
        end

    end
    
%% Compute OFM between ROOF and LAND    
elseif dem(x,y) == 999 && dem(xx,yy) < 999 || dem(x,y) < 999 && dem(xx,yy) == 999 %%|| abs(dem(x,y) - dem(xx,yy))/s_grid >= 1      
    if pos_oflow <= deps
        dqq = 0;
    else
        h_ofm = pos_oflow - deps;
        if slop_fic >= 0
            dqq = 0.6 * 0.5 * sqrt(2 * 9.81) * h_ofm^(3/2);                                  % Using equation: q = c & w * sqrt(2 * g) * z_ofm ^ 1.5; c = 0.6, w = 2; g = 9.81
        else
            dqq = - 0.6 * 0.5 * sqrt(2 * 9.81) * h_ofm^(3/2);
        end
    end   
end

%% Using implicit scheme and modifed velocity to avoid a negative water depth at the next time step
        % according to the water depth at previous time step
if TSP2D == 1
velocity = dqq / s_grid^2;
max_v = (pos_oflow - deps)/deltaT; 
if velocity ~= 0
    if abs(velocity) > max_v    
        if slop_fic >= 0
            velocity = max_v;
        else
            velocity = -max_v;
        end
    end
end

%% Using implicit scheme
elseif TSP2D == 2               
    velocity = dqq / s_grid^2;      
end
end

