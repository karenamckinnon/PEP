function yearByYearFigs(C, T95prime, GHCND_clustered, summerTime, yrs, summer_doy, ...
	lag_range, fig_folder, predictand)
 
nyrs = length(yrs); 
T95summer = reshape(T95prime(ismember(GHCND_clustered.doy, summer_doy)),[length(summer_doy) nyrs]);

% add stippling when value is above 80th percentile
C80 = prctile(C, 80, 2);

for ct = 1:nyrs
	disp([num2str(ct)])
	idx = find(year(summerTime) == yrs(ct));
	CUse = C(:, idx);
	
	yrIdx = find(year(GHCND_clustered.time) == yrs(ct));
	GHCND_clustered_short = GHCND_clustered;
	GHCND_clustered_short.time = GHCND_clustered_short.time(yrIdx);
	GHCND_clustered_short.doy = GHCND_clustered_short.doy(yrIdx);
	 
	predictandUse = predictand(year(predictand) == yrs(ct));
	%{
	if ~isempty(predictandUse)
		% calculate ROC score for each lead time
		[roc_score] = ...
			calcROC(GHCND_clustered_short, CUse, [], lag_range, predictandUse, ...
			summer_doy, 0, fig_folder, 'Pacific');
	else
		roc_score = NaN(size(lag_range));
	end
	%}
	 
	clf
	p1 = subplot(3,3,[1:6]);
	imagesc(summerTime(idx), lag_range, CUse);
	% add scatter where there are values > 80th percentile
	above80 = CUse > repmat(C80, [1 60]);
	[X Y] = meshgrid(summerTime(idx), lag_range);
	hold on
	plot(X(above80), Y(above80), '.', 'color', 0.5*[1 1 1])
	set(gca,'ydir','normal')
	ylim([0 50])
	xlim([summerTime(idx(1)) summerTime(idx(end))])
	set(gca,'xtick',[])
	colormap(flipud(lbmap(21,'redblue')))
	caxis([-2.5 2.5])
	set(gca,'fontsize',16)

	%{
	subplot(3,4,[4,8])
	plot([0.6 0.6],[min(lag_range) max(lag_range)], '--', 'color', 0.8*[1 1 1],'linewidth',2)
	hold on
	plot(roc_score, lag_range, 'k','linewidth',2)
	ylim([0 50])
	xlim([0 1])
	%}
	
	subplot(3,3,[7:9])
	T95plot = T95summer(:, ct);
	tPlot = summerTime(idx);
	cutoffHi = std(T95summer(:)) + mean(T95summer(:));

	% interpolate to fine resolution for plot
	X = linspace(tPlot(1),tPlot(end), 601);
	T95plotinterp = interp1(tPlot, T95plot, X);

	idxHi = T95plotinterp > cutoffHi;

	valsHi = T95plotinterp;
	valsHi(~idxHi) = cutoffHi;
	plot([tPlot(1) tPlot(end)],[cutoffHi cutoffHi],'--','color',0.8*[1 1 1],'linewidth',2);
	hold on
	h = patch([X fliplr(X)],[cutoffHi*ones(size(valsHi)) fliplr(valsHi)],'r');
	set(h,'edgecolor','none');
	hold on
	set(h,'edgecolor','none');
	plot(tPlot,T95plot,'-k','linewidth',2)
	
	xticks = tPlot(1):14:tPlot(end);
	set(gca,'xtick',xticks)
	set(gca,'xticklabel',xticks)
	datetick('x', 'mm/dd','keeplimits','keepticks')
	
	set(gca,'box','on')
	set(gca,'layer','top')
	set(gca,'fontsize',16)
	
	xlim([tPlot(1) tPlot(end)])
	ylim([-2 12])

	% calculate intraseasonal correlation
	% must be higher than 0.295 to be sig at the 95% level
	%rho = xcPH(CAvg, T95summer(:, ct), 1);
	%subplot(p1)
	%title(['Intra-seasonal correlation: ' num2str(1/100*round(rho*100)) ''])
	
	orient landscape

	figname = ['exampleFor' num2str(yrs(ct)) '.pdf'];
	% print('-dpng',[fig_folder '/' figname])

	set(gcf, 'Color', 'w')
	export_fig([fig_folder '/' figname])


end

% make a colorbar
clf
h = colorbar;
set(h,'fontsize',24)
colormap(flipud(lbmap(21,'redblue')))
caxis([-2.5 2.5])

set(gcf, 'Color', 'w')
export_fig([fig_folder '/colorbar.pdf'])