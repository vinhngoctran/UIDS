%============================================================================
%   coordinate.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program uses to convert UTM coordinates into Lat/Lon vectors (WGS84)
%============================================================================

function  [Lat,Lon] = coordinate(x,y,utmzone)

if (utmzone(1,4)>'X' || utmzone(1,4)<'C')
    fprintf('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
end
if (utmzone(1,4)>'M')
    hemis='N';
else
    hemis='S';
end
zone=str2double(utmzone(1,1:2));
sa = 6378137.000000 ; sb = 6356752.314245;
e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
e2cuadrada = e2 ^ 2;
c = ( sa ^ 2 ) / sb;
X = x - 500000;
if hemis == 'S' || hemis == 's'
    Y = y - 10000000;
else
    Y = y;
end
S = ( ( zone * 6 ) - 183 );
lat =  Y / ( 6366197.724 * 0.9996 );
v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
a = X / v;
a1 = sin( 2 * lat );
a2 = a1 * ( cos(lat) ) ^ 2;
j2 = lat + ( a1 / 2 );
j4 = ( ( 3 * j2 ) + a2 ) / 4;
j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
alfa = ( 3 / 4 ) * e2cuadrada;
beta = ( 5 / 3 ) * alfa ^ 2;
gama = ( 35 / 27 ) * alfa ^ 3;
Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
b = ( Y - Bm ) / v;
Epsi = ( ( e2cuadrada * a^ 2 ) / 2 ) * ( cos(lat) )^ 2;
Eps = a * ( 1 - ( Epsi / 3 ) );
nab = ( b * ( 1 - Epsi ) ) + lat;
senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
Delt = atan(senoheps / (cos(nab) ) );
TaO = atan(cos(Delt) * tan(nab));
longitude = (Delt *(180 / pi ) ) + S;
latitude = ( lat + ( 1 + e2cuadrada* (cos(lat)^ 2) - ( 3 / 2 ) * e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) ) * ( TaO - lat ) ) * (180 / pi);
Lat=latitude;
Lon=longitude;
end