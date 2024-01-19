%============================================================================
%   qmanhole.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate inflow to manhole of 1D SWMM5
%============================================================================

function [q_to_manhole of_rate] = qmanhole(deltaT,s_grid,deps,oflow,of_rate,node_depth_at_t,dem,node_invert)

if deltaT <= 1                                      %If deltaT < 1 second
    pos_oflow = oflow + of_rate * deltaT;
    
    if pos_oflow < 0
        pos_oflow = 0;
    elseif pos_oflow <= deps && of_rate < 0 && oflow >= deps
        pos_oflow = deps;
    end
    
    if pos_oflow <= deps
        q_to_manhole = 0;
        
    else
        z_node = node_depth_at_t + dem - node_invert;
        h_ofm = pos_oflow - deps;
        if z_node < dem
            q_to_manhole = 1 * 1.5 * sqrt(2 * 9.81) * h_ofm^(3/2);                                  % Using equation: q = c & w * sqrt(2 * g) * z_ofm ^ 1.5; c = 1, w = 1.5; g = 9.81
        elseif z_node >= dem && z_node <= h_ofm + dem
            q_to_manhole = 1 * 1.5 * 1.5 * sqrt(2 * 9.81) * (h_ofm + dem - z_node)^(1/2);          % Using equation: q = c & A * sqrt(2 * g) * sqrt(z_ofm + dem - z_node); A =  1.5 * 1.5
        else
            q_to_manhole = 0;
        end
        if q_to_manhole * deltaT / s_grid^2 >= h_ofm
            q_to_manhole = h_ofm * s_grid^2 / deltaT;
        end
    end
    of_rate = of_rate - q_to_manhole / s_grid^2;
    
else                                              % If no, we should normalize the outflow
    pre_oflow = oflow;
    for s_t = 1:deltaT
        pre_of_rate = of_rate;
        pos_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
        
        if pos_oflow < 0
            pos_oflow = 0;
        elseif pos_oflow <= deps && pre_of_rate < 0 && pre_oflow >= deps
            pos_oflow = deps;
        end
        
        if pos_oflow <= deps
            q_to_manhole = 0;
            
        else
            z_node = node_depth_at_t + dem - node_invert;
            h_ofm = pos_oflow - deps;
            if z_node < dem
                q_to_manhole = 1 * 1.5 * sqrt(2 * 9.81) * h_ofm^(3/2);                                  % Using equation: q = c & w * sqrt(2 * g) * z_ofm ^ 1.5; c = 0.6, w = 2; g = 9.81
            elseif z_node >= dem && z_node <= h_ofm + dem
                q_to_manhole = 1 * 1.5 * 1.5 * sqrt(2 * 9.81) * (h_ofm + dem - z_node)^(1/2);          % Using equation: q = c & A * sqrt(2 * g) * sqrt(z_ofm + dem - z_node);  A =  1.5 * 1.5;
            else
                q_to_manhole = 0;
            end
        end
        preq(s_t) = q_to_manhole;
        pre_of_rate = pre_of_rate - preq(s_t) / s_grid^2;
        pre_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
    end
    q_to_manhole = mean(preq);
    if q_to_manhole ~= 0
        if q_to_manhole * deltaT / s_grid^2 >= (oflow + of_rate * deltaT - deps)
            q_to_manhole = (oflow + of_rate * deltaT - deps) * s_grid^2 / deltaT;
        end
    end
    of_rate = of_rate - q_to_manhole / s_grid^2;
end