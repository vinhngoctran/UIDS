function meshgen()

DEM_DAT = textread(filename);
shp = DEM_DAT(1:n_row,1:n_col); 

%% Land: shpe == 0

%% Chanel: shp == 1

%% Oulet: shp == 2

%% Critical boundary: shp == 3

%% Soft boundary: shp == 4

%% Open boundary: shp == 5

%% Housing: shp == 6


%% Boundary: shp == 5
for x = 1:size(shp,1)-1
    for y = 1:size(shp,2)-1
        if shp(x,y) ~= -9999 &&  shp(x+1,y) == -9999 || shp(x,y) ~= -9999 &&  shp(x,y+1) == -9999 
            shp(x,y) = 5;
        elseif shp(x,y) == -9999 &&  shp(x+1,y) ~= -9999 
            shp(x+1,y) = 5;
        elseif shp(x,y) == -9999 &&  shp(x,y+1) ~= -9999
            
            shp(x,y+1) = 5;
        end
        
    end
end



end