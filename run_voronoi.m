function [stn_areas,stn_polylon,stn_polylat] = run_voronoi(varargin)
if ischar(varargin{1})
    demo = 1;
else
    demo = 0;
end

% Better version: simply computes the _full_ spherical voronoi patches and
% then gets the intersection via polybool. A bit slower, but vastly
% superior in terms of stability.
 
% Currently set up only for North America
% Use run_voronoi('demo') for a demo.

% Shapefile options: 
%  i) default matlab landareas -- excludes lots of islands!!!
%namer = shaperead('landareas','UseGeoCoords',true,'Selector',{@(name)strcmpi(name,'North and South America'),'Name'});
%  ii) manually constructed North America + surrounding islands, at
%  intermediate resolution:
%load shapefiles/CONUS_all.mat
load shapefiles/CONUS_smallerbuffer.mat
% 
% disp('Generating buffer zones...')
% latlim = [13 89];
% lonlim = [-179 -30];
% [lat_namer,lon_namer] = maptrimp(namer.Lat,namer.Lon,latlim,lonlim);
% [rpolylat,rpolylon,~,~] = reducem(lat_namer',lon_namer',0.1);
% disp('Precision buffer...')
% [lat_namer_b,lon_namer_b] = bufferm(rpolylat, rpolylon, 0.1, 'out'); % Buffer points for lat/lon precision of coastal stations.
% [rpolylon_b,rpolylat_b] = polybool('union',rpolylon,rpolylat,lon_namer_b,lat_namer_b);
% [lon_namer_b,lat_namer_b] = poly2cw(rpolylon_b,rpolylat_b); % Final clockwise-oriented mask to compare against. Can do similar for full world map.
% disp('Coarse bufer...')
% [latb,lonb] = bufferm(rpolylat, rpolylon, 30.0, 'out'); % Buffer points
% [rlatb,rlonb] = reducem(latb, lonb, 0.03);
% breakpoints = [0; find(isnan(rlatb))]; % Get largest closed polygon to close entire region.
% for j = 1:length(breakpoints)-1
%     segment_area(j) = areaint(rlatb(breakpoints(j)+1:breakpoints(j+1)-1),rlonb(breakpoints(j)+1:breakpoints(j+1)-1));
% end
% [~,idx] = sort(segment_area);
% rlatb = rlatb(breakpoints(idx(end))+1:breakpoints(idx(end)+1)-2); % (Cut off last point since it is identical to the first)
% rlonb = rlonb(breakpoints(idx(end))+1:breakpoints(idx(end)+1)-2);
%[rlonb,rlatb] = largest_poly(rlonb,rlatb);
%[rlatb,rlonb] = flatearthpoly(rlatb,rlonb);
[lat_namer_b,lon_namer_b,~,~] = reducem(namer.Lat',namer.Lon',0.1);
xyz_b = lonlat2ecef(rlonb,rlatb);
disp('Generating anchor points outside of buffer zone')
% Also generate quasi-uniform points outside of buffer zone:
xyz_out = randn(3,9000);
xyz_out = bsxfun(@rdivide, xyz_out, sqrt(sum(xyz_out.^2,1)));
[lon_out,lat_out] = ecef2lonlat(xyz_out);
IN = inpolygon(lon_out,lat_out,rlonb,rlatb);
lon_out = lon_out(~IN);
lat_out = lat_out(~IN);
xyz_out = lonlat2ecef(lon_out,lat_out);
if demo
    n = 15000;
    xyz = randn(3,n);
    xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));
    % Make sure points aren't too close
    [IDX,D] = knnsearch(xyz',xyz','K',2);
    [~,idx_remv] = find(D(:,2)<1e-03);
    xyz = xyz(:,setdiff(1:size(xyz,2),IDX(idx_remv,2)));
    [lon_rand,lat_rand] = ecef2lonlat(xyz);
    in = inpolygon(lon_rand,lat_rand,lon_namer_b,lat_namer_b);
    xyz_in = xyz(:,in);
    lon_in = lon_rand(in);
    lat_in = lat_rand(in);
else
    lon_in = varargin{1};
    lat_in = varargin{2};
    xyz_in = lonlat2ecef(lon_in,lat_in);
end
stn_polylon = cell(length(lon_in),1);
stn_polylat = cell(length(lon_in),1);
stn_areas = NaN(size(lon_in));
disp('Computing convex hull')
[P, K] = voronoisphere_2([xyz_in xyz_out],'resolution',0.0005*pi/180.0);
%[P, K] = voronoisphere_2([xyz_in xyz_b],'resolution',0.001*pi/180.0);
%[P, K] = voronoisphere_2(xyz_in,'resolution',0.001*pi/180.0);
[P_lon,P_lat] = ecef2lonlat(P);

if demo
figure(1)
clf
end
% Get only intersection with land areas.
disp('Finding intersection of polygons with land areas')
for j = 1:size(xyz_in,2)
    [lon_cw,lat_cw] = poly2cw(P_lon(K{j}),P_lat(K{j}));
    [stn_polylon{j},stn_polylat{j}] = polybool('intersection',lon_namer_b,lat_namer_b,lon_cw,lat_cw);
    try
    stn_areas(j) = sum(areaint(stn_polylat{j},stn_polylon{j})); % Units: fraction of surface of sphere
    catch err
        stn_areas(j) = NaN;
        j
        disp('Station not within CONUS');
    end
    if demo
        [LAT_union,LON_union] = polysplit(stn_polylat{j},stn_polylon{j});
        plot(lon_namer_b,lat_namer_b,'r-')
        hold on
        plot(rlonb,rlatb,'k-')
        for k = 1:length(LAT_union)
            patch(LON_union{k},LAT_union{k},[.6 .6 .6])
        end
        plot(lon_in(j),lat_in(j),'b.','MarkerSize',12)
        plot(lon_cw,lat_cw,'r.-','LineWidth',2)
        hold off
        title(['Area = ' num2str(stn_areas(j)*100,2) '% of globe']);
        %pause
    end
end 

if demo
   keyboard
end

return