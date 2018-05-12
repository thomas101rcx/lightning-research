%%%
function Dists = spheric_distance(lat1, lon1, lat2, lon2)
%%%
%%% This function calculate the spheric distance between two location on
%%% earth.   kilometers
%%%

format long;
lat1 = double(lat1);
lon1 = double(lon1);
lat2 = double(lat2);
lon2 = double(lon2);

R = 6371.0;
pdeg2rad = pi/180.0;

lat1 = lat1*pdeg2rad;
lat2 = lat2*pdeg2rad;
lon1 = lon1*pdeg2rad;
lon2 = lon2*pdeg2rad;

Dists = R*acos(sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(lon2-lon1));

