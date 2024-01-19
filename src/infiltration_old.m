%============================================================================
%   infiltration.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate infiltration rate by Green-Amplt method
%============================================================================

% General variables
% Variables for Green-Amplt method
%pa_infil: ARRAY PROVIDES INFORMATION ABOUT SOIL PARAMETERS
%inf_rate: RATE OF INFITRATION
%soil: TYPE OF SOIL
%d_oflow: EXCESS RAINFALL IN DEPTH UNIT(MM/deltaT,DEPTH OF EXCESS WATER LATER)
%deltaT: TIME INTERVAL
%du_theta: soil moisture deficit in the upper soil recovery zone
%d_theta: soil moisture deficit at the start of the current rainfall event

%Ks: HYDRAULIC CONDUCTIVITY [m/s]
%Wh: Suction head at the wetting front [m]
%SMD: MAXIMUM SOIL MOISTURE DEFICIT
%LU: Depth of upper soil recovery zone [m]
%Kr: Moisture deficit recovery constant [1/hour]

%%
function [d_theta, du_theta, v_infil, of_rate]=infiltration_old(x,y,deltaT,soil,d_theta,du_theta,v_infil,pa_infil,oflow,of_rate,imp,varargin)

if soil(x,y) ~= -9999
    iinf = soil(x,y);
else
    iinf=1;
end

% Read infitration parameters
Ks = pa_infil(iinf,1);              
Wh = pa_infil(iinf,2);
smd = pa_infil(iinf,3);
LU = 4 * sqrt(Ks*3600);
Kr = sqrt(Ks*3600) / 75;

if imp(x,y) < 0
    imp(x,y)=0.;
end

% Compute the infiltrated rainfall based on impervious area ratio
id_oflow = of_rate * (1 - imp(x,y)) + oflow/deltaT;

% If the available rainfall rate is zero (ia = 0) then:
if id_oflow == 0
    inf_rate = 0;
    delta_theta = Kr * smd * deltaT;
    du_theta = du_theta + delta_theta;
    v_infil = v_infil - delta_theta * LU;
    id_oflow = 0;
    
elseif id_oflow <= Ks    % If the available rainfall rate does not exceed the saturated hydraulic conductivity (ia ≤ Ks) then:
    inf_rate = id_oflow;
    v_infil = v_infil + inf_rate * deltaT;
    du_theta = du_theta - inf_rate * deltaT / LU;
    id_oflow = 0;
    
else  %If the available rainfall rate id_oflow exceeds Ks then: 
    Full_v = (Ks * Wh * d_theta) / (id_oflow - Ks);
    
    if v_infil >= Full_v
        v_1 = v_infil;
%         kc = Ks * deltaT + v_1 - Wh * d_theta * log(v_1 + Wh * d_theta);
        v_2 = v_1 + Ks * deltaT;
%         v_2 = fzero(@(v_2)kc +  Wh * d_theta * log(v_2 + Wh * d_theta),v_1);
        d_v = v_2 - v_infil;
        if d_v >= id_oflow * deltaT
            v_infil = v_infil + d_v;
            du_theta = du_theta - d_v / deltaT;
            inf_rate = id_oflow;
            id_oflow = 0;
        else
            d_v = v_2 - v_infil;
            v_infil = v_2;
            du_theta = du_theta - d_v / LU;
            inf_rate = d_v / deltaT;
            id_oflow = id_oflow - inf_rate;
        end
        
    elseif  v_infil + id_oflow * deltaT < Full_v
        inf_rate = id_oflow;
        v_infil = v_infil + inf_rate * deltaT;
        du_theta = du_theta - inf_rate * deltaT / LU;
        id_oflow = 0;
    else
        v_1 = Full_v;
%         kc = Ks * deltaT + v_1 - Wh * d_theta * log(v_1 + Wh * d_theta);
        v_2 = v_1 + Ks * deltaT;             
        d_v = v_2 - v_infil;
        v_infil = v_2;
        du_theta = du_theta - d_v / LU;
        inf_rate = d_v / deltaT;
        id_oflow = id_oflow - inf_rate;
    end
end
    
of_rate = of_rate - inf_rate ;

% if du_theta < 0 -> du_theta = 0  and   du_theta > SMD  -> du_theta = SMD
if du_theta < 0
    du_theta = 0;
elseif du_theta > smd
    du_theta = smd;
end  
end 

