%============================================================================
%   visualize.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is to visualize the model results
%============================================================================

function visualize(Exp_flo,head_plot,t,deltaT,dem,n_row,n_col)

for x = 1:n_row
    for y = 1:n_col
        if dem(x,y) == -9999 || dem(x,y) == 999
            Exp_flo(x,y) = NaN;
        end
    end
end

figure1=figure;set(gcf, 'Position', [300, 100, 600, 700]);

Time = t*deltaT/3600; Name = ['Inundation map at t = ',num2str(round(Time,1)),' hour'];

title(Name);hold on

mapshow(Exp_flo,head_plot,'DisplayType','surface')

xlabel('x (easting in meters)');

ylabel('y (northing in meters)');

demcmap(Exp_flo);

caxis([0 1]);

colorbar;

end