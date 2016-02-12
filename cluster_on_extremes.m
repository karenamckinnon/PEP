function GHCND_clustered = cluster_on_extremes(GHCND_QC, maxclust, summer_doy, makeplots, ...
	station_latlims, station_lonlims, fig_folder, cacheDir)
% Cluster station data based on whether days are extreme on the same day
% 
% Input: GHCND_QC = structure containing info about stations that passed QC
%		 note that anom contains all data, not just summer
%        maxclust = number of clusters to use. Note that some stations end
%        up being grouped into their own cluster, so this is larger than
%        the actual numer of clusters being used in the end
%		 summer_doy = days to use for summer
% Output: GHCND_clustered = structure amended to include the cluster id


% anomalies from baseline during summer
summer_anom = GHCND_QC.anom(:, ismember(GHCND_QC.doy, summer_doy));

% Calculate binary field
binary_temperature = summer_anom;
binary_temperature(binary_temperature <= 0) = 0;
binary_temperature(binary_temperature > 0) = 1;

% cluster on above or below baseline
rng('default');
cluster_input.binary_temperature = binary_temperature;
cluster_input.maxclust = maxclust;
cluster_input.dist_method = 'jaccard';
cluster_input.linkage_method = 'average';

% Hash to see if this has already been done...
cluster_hash = DataHash(cluster_input);
plot_dendro = 0;
plot_silhouette = 0;

figpath = fig_folder;
if exist([cacheDir '/cluster' cluster_hash '.mat'],'file') % load existing?
    disp('Clustering already computed... loading')
    load([cacheDir '/cluster' cluster_hash '.mat'])
else % If not, cluster and save:
    disp('First time clustering this data... will take a while.')
    tic
    [clusters,c,Z] = cluster_qdata(binary_temperature, cluster_input.maxclust, ...
    	cluster_input.dist_method, cluster_input.linkage_method, plot_dendro, ...
    	plot_silhouette, figpath);
    toc
    timestamp = datestr(now);
    save([cacheDir '/cluster' cluster_hash '.mat'],'clusters','c','Z','timestamp')
end


% make plot of different numbers of clusters?

clusterInd = maxclust;

if makeplots
	for ii = 1:maxclust
		cluster_temp = clusters(:, ii);
		cluster_temp = relabel_clusters(cluster_temp);
		map_clustered_stations(GHCND_QC, cluster_temp, fig_folder, ...
			maxclust, station_latlims, station_lonlims)
	end
end

clusters = clusters(:,clusterInd);
[clusters,cluster_counts] = relabel_clusters(clusters);

GHCND_clustered = GHCND_QC;
GHCND_clustered.clusterid = clusters;

return




