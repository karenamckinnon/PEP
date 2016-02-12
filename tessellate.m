function [station_areas, stn_polylon, stn_polylat] = tessellate(GHCND_clustered, cacheDir)
% Calculate voronoi tessellation for station locations
% Result used to estimate areal averages based on point data at stations
%
% Input: GHCND_clustered = structure containing station data
% Output: station_areas = area of each tessellation in degrees^2 (?)
%         stn_polylon = longitude of station (?)
%         stn_polylat = latitude of station (?)

rng('default');
voronoi_input.lats = GHCND_clustered.location(:, 2) + ...
    0.001*randn(size(GHCND_clustered.location(:, 2)));
voronoi_input.lons = GHCND_clustered.location(:, 1);
 
area_hash = DataHash(voronoi_input);
if exist([cacheDir '/voronoi' area_hash '.mat'],'file') % load existing?
    disp('Voronoi tessellation already computed... loading')
    load([cacheDir '/voronoi' area_hash '.mat'])
else % If not, cluster and save:
    disp('First time computing tessellation on this list... ')
    [areas,polylon,polylat] = run_voronoi(voronoi_input.lons,voronoi_input.lats);
    timestamp = datestr(now);
    save([cacheDir '/voronoi' area_hash '.mat'], ...
    	'areas', 'polylon', 'polylat', 'timestamp');
end
station_areas = areas;
stn_polylon = polylon;
stn_polylat = polylat;

return