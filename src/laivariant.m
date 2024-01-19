%============================================================================
%   cflow.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program calculates the daily or sub-daily value of LAI sing a continuous, 
%   piecewise linear function of LAI_max and LAI_min.
% 
%   Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998),Crop 
%   evapotranspiration:  Guidelines for computing crop requirements, FAO 
%   Irrigation and Drainage Paper 56, Food and Agriculture Organization 
%   (FAO), Rome, Italy
%============================================================================

function [LAI_d,LAI_m] = laivariant(LAI_max,LAI_min,time)

% INPUTS:
%       LAI_min         Reduction in leaf area index in the dormant season
%                           compared to the growing season [%]                    
%       LAI_max         Annual maximum leaf area index [m2 m-2]              
%                           time time vector of date numbers of the period for  
%                           which leaf area index will be calculated          
%                           (dates at the beginning of the time steps)
% 
% OUTPUTS:
%       LAI_d           daily/sub-daily LAI [m2 m-2] 
%       LAI_m           monthly LAI [m2 m-2]                               
% 
% 
% NOTE: this function assumes that LAI is equal to its maximum value 
% in June, July and August, and to its minimum value in December, January
% and February.
% 
%--------------------------------------------------------------------------
% 1. Define the different seasons 
%--------------------------------------------------------------------------
winter = [12 1 2];      % dormant season (LAI is equal to its maximum value)
summer = [6 7 8];       % growing season (LAI si equal to its minimum value)
spring = [3 4 5]; 
autumn = [9 10 11];

%--------------------------------------------------------------------------
% 2. Assess monthly LAI
%--------------------------------------------------------------------------
LAI_m = nan(12,1);      % Initialise variable
LAI_m(winter) = LAI_min/100*LAI_max;
LAI_m(summer) = LAI_max;
LAI_m(spring) = LAI_min/100*LAI_max/4*(6-spring)+LAI_max/4*(spring-2);
LAI_m(autumn) = LAI_min/100*LAI_max/4*(autumn-8)+LAI_max/4*(12-autumn);

%--------------------------------------------------------------------------
% 3. Assess LAI
%--------------------------------------------------------------------------
LAI_d = nan(length(time),1); % Initialise variable
date_vec = datevec(time);

for m=1:12
     idx_m = date_vec(:,2)==m;
     LAI_d(idx_m) = LAI_m(m);
end
end
