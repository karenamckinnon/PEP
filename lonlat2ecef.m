% LONLAT2ECEF - convert latitude, longitude to
%            cartesian coordinates using the wgs84 ellipsoid but on the
%            unit sphere.
% 
% USAGE:
% xyz = lonlat2ecef(lat,lon)
% 
% x = ECEF X-coordinate
% y = ECEF Y-coordinate
% z = ECEF Z-coordinate
% lat = geodetic latitude (degrees) on [-90,90]
% lon = longitude (degrees) on [-180 180]

function xyz=lonlat2ecef(lon,lat)
lon(lon<0) = lon(lon<0) + 360;
% WGS84 ellipsoid constants:
a = 1;% Unit sphere... for earth use 6378137;
e = 8.1819190842622e-2; 
N = a ./ sqrt(1 - e^2 .* sin(lat*pi/180.0).^2);

% results:
xyz = NaN(3,length(lon));
xyz(1,:) = (N) .* cos(lat*pi/180.0) .* cos(lon*pi/180.0);
xyz(2,:) = (N) .* cos(lat*pi/180.0) .* sin(lon*pi/180.0);
xyz(3,:) = ((1-e^2) .* N) .* sin(lat*pi/180.0);

return