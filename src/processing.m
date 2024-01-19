%============================================================================
%   processing.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is to process the data according to selected deltaT
%============================================================================

function data = processing(dataraw, deltaT,dt,t_dur,type)

% Check the size of raw data
[x,y,z] = size(dataraw);

% Generate time series data
k = 0;
if z == 1
    for i = 1:t_dur/dt
        for j = 1:dt/deltaT
            k = k + 1;
            if type == 1                    % For cummulative precipitation
                data(k,:) = dataraw(i,:)/(dt/deltaT);
            elseif type == 0                % Constant
                data(k,:) = dataraw(i,:);
            elseif type == 2                % Linear interpolation
                data(k,:) = dataraw(i,:)+(dataraw(i+1,:)-dataraw(i,:))*j/(dt/deltaT);
            end
            
        end
    end
    
% Generate gridded data
else
    for i = 1:t_dur/dt
        for j = 1:dt/deltaT
            k = k + 1;
            if type == 1                    % For precipitation
                data(:,:,k) = dataraw(:,:,i)/(dt/deltaT);
            elseif type == 0                % Constant
                data(:,:,k) = dataraw(:,:,i);
            elseif type == 2                % Linear interpolation
                data(:,:,k) = dataraw(:,:,i)+(dataraw(:,:,i+1)-dataraw(:,:,i))*j/(dt/deltaT);
            end
        end
    end
end
end