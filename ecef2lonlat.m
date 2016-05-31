% ECEF2LONLAT - convert earth-centered earth-fixed (ECEF)
%            cartesian coordinates to latitude, longitude,
%
% USAGE:
% [lon,lat] = ecef2lonlat(xyz)
%
% lat = geodetic latitude (degrees)
% lon = longitude (degrees)
% x = ECEF X-coordinate 
% y = ECEF Y-coordinate 
% z = ECEF Z-coordinate 

function [lon,lat] = ecef2lonlat(xyz)
x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);

a = 1; % Unit sphere... use 6378137 if on earth...;
e = 8.1819190842622e-2;

b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);
lon = 180.0/pi*atan2(y,x);
lat = 180.0/pi*atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));

lon = mod(lon,360.0);
lon(lon>180.0) = lon(lon>180.0)-360;

return