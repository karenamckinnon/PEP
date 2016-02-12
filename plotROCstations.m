function roc_scoreStation = plotROCstations(GHCND_clustered, cutofflength, C, lag_range, lagsToPlot, ...
	fig_folder, summerDays, yrs, summer_doy, cone_ratio, predictandName)
 
% Station by station
rng('default') 
clear cache_input
% also calculate for each station
stationAnom = bsxfun(@minus, GHCND_clustered.anomTotal, nanmean(GHCND_clustered.anomTotal(:, summerDays), 2));
cache_input.data = stationAnom;
cache_input.lag_range = lag_range;
cache_input.predictand = predictandName;
ROCstationHash = DataHash(cache_input);
matname =['/n/huybers_lab/common/ghcnd/analysis/composites/cache/ROC_' ROCstationHash '.mat'];
if exist(matname, 'file')
	load(matname)
else
	for ct = 1:size(stationAnom, 1)
		disp(['' num2str(ct) ''])
		[startDateStation{ct}, endDateStation{ct}, totalTempStation{ct}, alldaysStation{ct}] = ...
			findHotDays(GHCND_clustered, cutofflength, stationAnom(ct, :), summerDays, ...
			yrs, summer_doy);
	end

	roc_scoreStation = NaN(length(lag_range), size(stationAnom, 1));
	predictand = eval(cache_input.predictand);
	for ct = 1:size(stationAnom, 1)
		disp(['' num2str(ct) ''])
		[roc_scoreStation(:, ct)] = ...
			calcROC(GHCND_clustered, C, cone_ratio, lag_range, predictand{ct}, ...
				summer_doy, 0, fig_folder, GHCND_clustered.id(ct, :));


	end
	timestamp = datestr(now);
	save(matname, 'roc_scoreStation', 'timestamp');
end
 
cmap = varycolor(5);
cmap = [0.3*[1 1 1];cmap];
for jj = 1:length(lagsToPlot)
	idx = find(lag_range == lagsToPlot(jj));
	clf
	m_proj('miller', 'lat', [25 50], 'lon', [-105 -60]);
	hold on
	m_coast('patch', 0.9*[1 1 1],'edgecolor','none');
	m_grid('linestyle','none','fontsize',24,'tickdir','in');
	[sx, sy] = m_ll2xy(GHCND_clustered.location(:, 1), GHCND_clustered.location(:, 2));
	colormap(cmap)
	caxis([0.45 0.75])
	h = colorbar;
	set(h,'fontsize',24)
	hold on
	scatter(sx, sy, 150, roc_scoreStation(idx,:), '.','linewidth',2)
	title(['Lead time: ' num2str(lagsToPlot(jj)) ' days'])
	set(gca,'fontsize',14)
	orient landscape
	set(gcf,'Renderer','zbuffer')
	fig_name = ['ROCScores_StationMaps_lag' num2str(lagsToPlot(jj)) '_' predictandName '.ps'];
	print('-dpsc',[fig_folder '/' fig_name]);
end

% Make histograms
xi = linspace(0.3,0.8,100);
bw = 0.02;
clf
hold on
colors = flipud(lbmap(length(lagsToPlot),'brownblue'));
for jj = 1:length(lagsToPlot)
	idx = find(lag_range == lagsToPlot(jj));
	fi = ksdensity(roc_scoreStation(idx,:),xi,'bandwidth',bw);
	h(jj) = plot(xi,fi,'color',colors(jj,:),'linewidth',4);
	legend_text{jj} = [num2str(lagsToPlot(jj)) ' days'];
end
grid on
hl = legend(h,strtrim(cellstr(legend_text)),'location','northwest');
set(hl,'fontsize',20)
xlabel('ROC score')
ylabel('Probability density')
xlim([min(xi) max(xi)])
figname = ['ROCScores_StationDensity_' predictandName '.ps'];

orient landscape
print('-dpsc',[fig_folder '/' figname])

