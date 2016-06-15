function map_clustered_stations(GHCND_QC, cluster_temp, fig_folder, ...
	maxclust, station_latlims, station_lonlims)
 
clf
m_proj('miller','lat',station_latlims,'lon',station_lonlims);
[sx, sy] = m_ll2xy(GHCND_QC.location(:, 1), GHCND_QC.location(:, 2));
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy');
% try plotting all
plotMax = 10;
idx = cluster_temp <= plotMax; % only plot the five largest clusters
cluster_tempPlot = cluster_temp(idx);
sx = sx(idx);
sy = sy(idx);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy')
C = [rgb('cornflowerblue');...
	rgb('navajowhite');...
	rgb('darkorchid');...
	rgb('deeppink');... 
	rgb('lightseagreen');...
	rgb('teal');...
	rgb('firebrick');...
	rgb('darkslateblue');...
	rgb('gold');...
	rgb('black')];
if size(C,1) < max(cluster_tempPlot),
	C = varycolor(max(cluster_tempPlot));
end
C = C(1:max(cluster_tempPlot),:);
colormap(C)
hold on
scatter(sx, sy, 150, cluster_tempPlot, '.','linewidth',2)
set(gcf,'Renderer','zbuffer')
orient landscape
print('-dpsc',[fig_folder '/cluster_map_for_' num2str(max(cluster_temp)) '.ps'])