%============================================================================
%   runoffsimp.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is to simulate runoff with a conceptual
%   mesoscale Hydrologic Model (mHm)
%============================================================================



function [runoff,baseflow,slowflow,fastflow,rstorage,lstorage,perc] = runoffsimp...
                    (pa_runoff,runof_reg,infilrate,rstorage, lstorage,surfrunoff)

%%
% INPUT
% infilrate:            Input to the soil layer [m d-1]
% k0                    Recession coefficient of the upper
%                               reservoir, upper outlet [d]
% k1                    Recession coefficient of the upper reservoir,
%                               lower outlet [d]
% k2                    Baseflow recession coefficient [d]
% kp                    Percolation coefficient [d]
% alpha                 Exponent for the upper reservoir [-]
% um                    Threshold water depth in upper reservoir
%                               (for Runoff contribution) [m]
% karst_loss            Karstic percolation loss [-]
%     
% OUTPUT
% lstorage              Upper soil storage [m]
% rstorage              Groundwater storage [m]
% perc                  Percolation [m s-1]
% fastflow              Fast runoff component [m s-1]
% slowflow              Slow runoff component [m s-1]
% baseflow              Baseflow [m s-1]

%% Parameterization
k0 = pa_runoff(runof_reg,1);
k1 = pa_runoff(runof_reg,2);
k2 = pa_runoff(runof_reg,3);
kp = pa_runoff(runof_reg,4);
alpha = pa_runoff(runof_reg,5);
um = pa_runoff(runof_reg,6);

% Runoff generation for the saturated zone
[lstorage, rstorage,slowflow,fastflow,perc] = unsatzone(k1, kp, k0, alpha,...
    infilrate, um,rstorage, lstorage);

% Runoff generation for the saturated zone
[rstorage, baseflow] = satzone(k2, rstorage);

runoff = baseflow + slowflow + fastflow + surfrunoff;
end


%%
function [lstorage, rstorage,slowflow,fastflow, perc] = unsatzone(k1, kp, k0, alpha,...
    infilrate, um,rstorage, lstorage)

% SOIL LAYER BETWEEN UNSATURATED AND SATURATED ZONE
lstorage = lstorage + infilrate;

% AST INTERFLOW WITH THRESHOLD BEHAVIOUR
fastflow   = 0;
if (lstorage > um)
    fastflow   = min( (1/(k0 * 86400) * (lstorage - um)), (lstorage));
end
    lstorage = lstorage -  fastflow;
    
% SLOW PERMANENT INTERFLOW
    slowflow= 0;
    
    if lstorage > 0
        slowflow= min((1/(k1 * 86400) *(lstorage^(1 + alpha) )),(lstorage));
    end
        lstorage = lstorage - slowflow;
        
% PERCOLATION FROM SOIL LAYER TO THE SATURATED ZONE
        perc  = 1/(kp*86400) * lstorage;
        
%  Taking into account for the KARSTIC aquifers
%  karstic loss gain or loss if Karstic aquifer is present...
        if(lstorage > perc )
            lstorage = lstorage - perc;
            rstorage = rstorage + perc ;
        else
            rstorage = rstorage + lstorage;
            rstorage = 0;
        end       
end

%%
function [rstorage, baseflow] = satzone(k2, rstorage)

    if (rstorage > 0)
       baseflow   = 1/(k2*86400) * rstorage;
       rstorage = rstorage - baseflow;
    else
       baseflow   = 0;
       rstorage = 0;
    end
end