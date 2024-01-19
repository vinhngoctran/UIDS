%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  Overland Flow-Urban Flood Model (OFM-Urban)
%                         © Vinh Ngoc Tran, 2021
%                      HYDROLAB, University of Ulsan
%
% **OFM-Urban** is an open-source Matlab-based software package for the
% simulation of hydrological processes in complex watershed systems.
%
% Version: 1.0
% Modified date: 2021
% 
%
% Copyright (C) 2021 Free Software Foundation, Inc. <http://fsf.org/>
% Everyone is permitted to copy and distribute verbatim copies
% of this license document, but changing it is not allowed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -------------------------------------------------------------------------%
clear all; clc; close all;
addpath(genpath('src'));
addpath(genpath('Examples'));
addpath(genpath('Results'));

%% TEST 01: Flooding a disconnect water body
clear all; close all; clc
Filename = 'Test_01.inf';
addpath(genpath('src'));
Results = ofm_urban(Filename);
% Results = ofm_urban_nosorting(Filename);
load('Results\Test01_test.mat')

for i=1:1200
   R(i) =  Results.Exp_flo(7,43,i)+dem(7,43);
   R2(i) = Results.Exp_flo(7,63,i)+dem(7,63);
end
for i=1:120
    RR(i) = mean(R((i-1)*10+1:i*10));
    RR2(i) = mean(R2((i-1)*10+1:i*10));
end
figure;
plot(RR,'linewidth',2); hold on;
plot(RR2,'linewidth',2);




%% TEST 02:  Filling of floodplain depressions
clear all; close all;clc;
addpath(genpath('src'));
Filename = 'Test_02.inf';
Results = ofm_urban(Filename);

load('Results\Test02.mat','Results','dem','inun_max','oflow','boundary_dat')
for i=1:2820
    l=0;
    for j=1:4
        for k=1:4
            l=l+1;
            R(i,l) =  Results.Exp_flo(14+25*(k-1),14+25*(j-1),i)+dem(14+25*(k-1),14+25*(j-1));
        end
    end
end
% for i=1:120
%     RR(i) = mean(R((i-1)*10+1:i*10));
%     RR2(i) = mean(R2((i-1)*10+1:i*10));
% end

for i=1:16
    figure;
plot(R(:,i),'linewidth',2);
pause
end

%% TEST 03:  Momentum conservation over a small obstruction
% ===> DO NOT USE THIS TEST BECAUSE OFM-Urban IS NOT A FULLY CONVERSE MOMENTUM BASED MODEL
clear all; close all;clc;
addpath(genpath('src'));
Filename = 'Test_03.inf';
Results = ofm_urban(Filename);
load('Results\Test03.mat','Results','dem','inun_max')

for i=1:900
   R(i) =  Results.Exp_flo(10,30,i)+dem(10,32);
   R2(i) = Results.Exp_flo(10,51,i)+dem(10,51);
end
figure;
plot(R,'linewidth',2); hold on;
plot(R2,'linewidth',2);

% for i=1:120
%     RR(i) = mean(R((i-1)*10+1:i*10));
%     RR2(i) = mean(R2((i-1)*10+1:i*10));
% end
figure;
plot(R,'linewidth',2); hold on;



%% TEST 04:  Valley flooding
clear all; close all;clc;
addpath(genpath('src'));
Filename = 'Test_05.inf';
Results = ofm_urban(Filename);
load('Results\Test05.mat','Results','dem','inun_max','WaterVec')
x_cor=232324.844;
y_cor=829983.063;
n_col=255;
n_row=225;
s_grid = 50;
outpoint = [235200,832400;
236700,833800;
237800,835200;
239400,838000;
243300,840300;
235700,832500;
237700,835500];

% Compute location of outlet
x_node_c = x_cor;
y_node_c = y_cor + n_row * s_grid;
for i=1:n_row
    for j=1:n_col
        x = x_node_c + j * s_grid - s_grid/2;
        y = y_node_c - i * s_grid + s_grid/2;  
        for k=1:size(outpoint,1)
            if abs(x-outpoint(k,1)) <= s_grid/2 && abs(y-outpoint(k,2)) <= s_grid/2
                x_out(k) = i;
                y_out(k) = j;
            end
        end
    end
end

for jj=1:7
for i=1:10800
   R(i,jj) =  Results.Exp_flo(x_out(jj),y_out(jj),i)+dem(x_out(jj),y_out(jj));
%    V(i,jj) =  WaterVec(x_out(jj),y_out(jj),i);
end
end

for i=1:7
figure;
plot(R(:,i),'linewidth',2); hold on;
pause
end

%% TEST 08A:  Rainfall and point source surface flow in urban areas
clear all; close all;clc;
addpath(genpath('src'));
Filename = 'Test_08A.inf';
Results = ofm_urban(Filename);

load('Results\Test08A.mat','Results','dem','inun_max')

% OUtput
x_cor=263976.000;
y_cor=664407.500;
n_col=483;
n_row=201;
s_grid = 2;
outpoint = [264682,664581;
264538,664665;
264356,664487;
264202,664553;
264334,664561;
264574,664555;
264710,664699;
264308,664647;
264222,664611];

% Compute location of outlet
x_node_c = x_cor;
y_node_c = y_cor + n_row * s_grid;
for i=1:n_row
    for j=1:n_col
        x = x_node_c + j * s_grid - s_grid/2;
        y = y_node_c - i * s_grid + s_grid/2;  
        for k=1:size(outpoint,1)
            if abs(x-outpoint(k,1)) <= s_grid/2 && abs(y-outpoint(k,2)) <= s_grid/2
                x_out(k) = i;
                y_out(k) = j;
            end
        end
    end
end

for jj=1:9
for i=1:300
   R(i,jj) =  Results.Exp_flo(x_out(jj),y_out(jj),i)+dem(x_out(jj),y_out(jj));
end
end
for i=1:9
figure;
plot(R(:,i),'linewidth',2); pause
end
%% TEST 08B:  Surface flow from a surcharging sewer in urban areas
clear all; close all;clc;
addpath(genpath('src'));
Filename = 'Test_08B.inf';
Results = ofm_urban(Filename);


load('Results\Test08B.mat','Results','dem','inun_max')
% OUtput
x_cor=263976.000;
y_cor=664407.500;
n_col=483;
n_row=201;
s_grid = 2;
outpoint = [264682,664581;
264538,664665;
264356,664487;
264202,664553;
264334,664561;
264574,664555;
264710,664699;
264308,664647;
264222,664611];

% Compute location of outlet
x_node_c = x_cor;
y_node_c = y_cor + n_row * s_grid;
for i=1:n_row
    for j=1:n_col
        x = x_node_c + j * s_grid - s_grid/2;
        y = y_node_c - i * s_grid + s_grid/2;  
        for k=1:size(outpoint,1)
            if abs(x-outpoint(k,1)) <= s_grid/2 && abs(y-outpoint(k,2)) <= s_grid/2
                x_out(k) = i;
                y_out(k) = j;
            end
        end
    end
end

for jj=1:9
for i=1:300
   R(i,jj) =  Results.Exp_flo(x_out(jj),y_out(jj),i)+dem(x_out(jj),y_out(jj));
end
end
for i=1:9
figure;
plot(R(:,i),'linewidth',2); pause
end

%% Application: Gangnam Area

% Gangnam 2010
Filename = 'GN2010.inf';
ofm_urban(Filename);

% Gangnam 2011
Filename = 'GN2011.inf';
ofm_urban(Filename);

% Gangnam 2013
Filename = 'GN2013.inf';
ofm_urban(Filename);

for t = 1:size(Res_flooddepth,3)
    for x = 1:n_row
        for y = 1:n_col
            if dem(x,y) == -9999 || dem(x,y) == 999
                Res_flooddepth(x,y,t) = NaN;
            end
        end
    end
end
figure1=figure;set(gcf, 'Position', [300, 100, 900, 1000]);
for t = 1:size(Res_flooddepth,3)
    HH = t*deltaT;
    Name = ['Inundation map at t = ',num2str(round(HH,1)),' second'];
    title(Name);hold on
    mapshow(Res_flooddepth(:,:,t),head_plot,'DisplayType','surface')
    xlabel('x (easting in meters)')
    ylabel('y (northing in meters)')
    demcmap(Res_flooddepth(:,:,t))
    caxis([0 0.3])
    cmp = colormap(parula);
    cmp = flipud(cmp);
    colormap(cmp);
    colorbar
    pause(0.01)   
end
