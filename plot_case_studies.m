function plot_case_studies(C, GHCND_clustered, ...
	sst_comp, hgt300_comp, lags_to_plot, lag_range, summer_doy, ...
	domainlatlims, domainlonlims, date_to_plot, fig_folder)
  
% make plots
station_lons =  GHCND_clustered.location(:, 1);
station_lats = GHCND_clustered.location(:, 2);
 
summer_idx = find(ismember(GHCND_clustered.doy, summer_doy));

% summer day index associated with case study
for ii = 1:length(date_to_plot)
	plot_idx = find(GHCND_clustered.time == date_to_plot(ii));
	sst_idx = find(sst_comp.time == date_to_plot(ii));
	z300_idx = find(hgt300_comp.time == date_to_plot(ii));
	for ct = 1:length(lags_to_plot);
		lag = lags_to_plot(ct);

		clf
		orient landscape
		set(gcf,  'PaperPositionMode', 'manual');
		set(gcf,  'PaperUnits', 'inches');
		set(gcf,  'PaperPosition', [0.25 0.25 8.5 11]);
		set(gcf,'visible','off')

		fig_name = ['case_study_for_' datestr(date_to_plot(ii)) '_lead_' num2str(lag) '_days.ps'];
		sst_mat = sst_comp.anom(:, :, sst_idx - lag);
		hgt300_mat = hgt300_comp.anom(:, :, z300_idx - lag);
		m_proj('miller','lat',[-10 80],...
			'lon',[120 330]);
		m_pcolorKM(sst_comp.lon, sst_comp.lat, sst_mat');
		shading interp
		m_coast('patch',[.7 .7 .7],'edgecolor','none');
		m_grid('box','fancy');
		hold on
		dk_turq = rgb('darkturquoise');
		% m_plot(station_lons+360,station_lats,'.','color',dk_turq)
		dk_green = rgb('darkgreen');
		m_plot(domainlonlims,[min(domainlatlims) min(domainlatlims)],'color',dk_green,'linewidth',2)
		m_plot(domainlonlims,[max(domainlatlims) max(domainlatlims)],'color',dk_green,'linewidth',2)
		m_plot([min(domainlonlims) min(domainlonlims)],domainlatlims,'color',dk_green,'linewidth',2)
		m_plot([max(domainlonlims) max(domainlonlims)],domainlatlims,'color',dk_green,'linewidth',2)
		cmap = flipud(lbmap(17, 'redblue'));
		cmap(9,:) = [1 1 1];
		colormap(cmap)
		h = colorbar;
		ylabel(h, 'SST anomaly (\circC)', 'fontsize', 14)
		clines = linspace(-400, 400, 17);
		m_contour(hgt300_comp.lon,hgt300_comp.lat,hgt300_mat',clines(clines>0),'linewidth',2,'col','k','linestyle','-');
		m_contour(hgt300_comp.lon,hgt300_comp.lat,hgt300_mat',clines(clines<0),'linewidth',2,'col','k','linestyle','--');
		caxis([-2 2])
		if ct == 1
			title(['' num2str(lag) ' day before ' datestr(date_to_plot(ii)) ''])
		else
			title(['' num2str(lag) ' days before ' datestr(date_to_plot(ii)) ''])
		end
		set(gcf,'Renderer','zbuffer')
		orient landscape
		print('-dpsc',[fig_folder '/' fig_name])
	end

end






