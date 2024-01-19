%============================================================================
%   qnode.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimates free outflow to sewer
%============================================================================

function q_to_freeoutlet = qnode(deltaT,of_rate,s_grid,deps,oflow)

if deltaT <= 1                                                  %If deltaT < 1 second
    pos_oflow = oflow + of_rate * deltaT;
    
    if pos_oflow <= deps
        q_to_freeoutlet = 0;
        
    else
        h_ofm = pos_oflow - deps;
        q_to_freeoutlet = 1 * s_grid * sqrt(2 * 9.81) * h_ofm^1.5;                      % Using equation: q = c & w * sqrt(2 * g) * z_ofm ^ 1.5; c = 1, w = 1; g = 9.81
        
        if q_to_freeoutlet * deltaT  / s_grid^2 >= h_ofm
            q_to_freeoutlet = h_ofm * s_grid^2 / deltaT;
        end
    end
    
else                                                            % If no, we should normalize the outflow
    pre_oflow = oflow;
    for s_t = 1:deltaT
        pre_of_rate = of_rate;
        pos_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
        
        if pos_oflow <= deps
            q_to_freeoutlet = 0;
            
        else
            h_ofm = pos_oflow - deps;
            q_to_freeoutlet = 1 * s_grid * sqrt(2 * 9.81) * h_ofm^1.5;                      % Using equation: q = c & w * sqrt(2 * g) * z_ofm ^ 1.5; c = 1, w = 1; g = 9.81
            
            if q_to_freeoutlet  / s_grid^2 >= h_ofm
                q_to_freeoutlet = h_ofm * s_grid^2;
            end
        end
        preq(s_t) = q_to_freeoutlet;
        pre_of_rate = pre_of_rate - preq(s_t) / s_grid^2;
        pre_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
    end
    q_to_freeoutlet = mean(preq);
    
    if q_to_freeoutlet ~= 0
        if q_to_freeoutlet * deltaT  / s_grid^2 >= (oflow + of_rate * deltaT - deps)
            q_to_freeoutlet = (oflow + of_rate * deltaT - deps) * s_grid^2 / deltaT;
        end
    end
end
end