%============================================================================
%   cflow.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   % Program is to compute discharge at the outlet
%============================================================================

function q_out = qoutlet(s_grid,deltaT,oflow,x_out,y_out,iout,of_rate,D_out,W_out,pflow)

if deltaT <= 1
    if oflow(x_out(iout),y_out(iout)) + of_rate(x_out(iout),y_out(iout)) * deltaT > D_out(iout)
        q_out = W_out(iout) * pflow * ((oflow(x_out(iout),y_out(iout)) + of_rate(x_out(iout),y_out(iout)) * deltaT - D_out(iout))^(5/3));
    else
        q_out = 0;
    end
else
    pre_oflow = oflow(x_out(iout),y_out(iout));
    for s_t = 1:deltaT
        pre_of_rate = of_rate(x_out(iout),y_out(iout));
        pos_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
        
        if pos_oflow <= D_out(iout)
            q_out = 0;
        else
            q_out = W_out(iout) * pflow * (pos_oflow - D_out(iout))^(5/3);
        end
        preq(s_t) = q_out;
        pre_of_rate = pre_of_rate - preq(s_t) / s_grid^2;
        pre_oflow = pre_oflow + pre_of_rate * deltaT/deltaT;
    end
    q_out = mean(preq);
    %	Limit the outflow by availability
    if  q_out ~= 0
        if q_out * deltaT / s_grid^2 >= (oflow(x_out(iout),y_out(iout))+ of_rate(x_out(iout),y_out(iout)) * deltaT - D_out(iout))
            q_out = (oflow(x_out(iout),y_out(iout))+ of_rate(x_out(iout),y_out(iout)) * deltaT - D_out(iout)) * s_grid^2 / deltaT;
        end
    end
end
end