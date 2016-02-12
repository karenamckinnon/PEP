function makeMaps(lag_range, lagsToPlotStart, lagsToPlotEnd, sst_comp, hgt300_comp, ...
	hgt850_comp, uwnd10m_comp, vwnd10m_comp, ...
	qnet_comp, alldays, WAF, fig_folder, ...
	domainlatlims, domainlonlims, sstCompSurrBnds, plot1, plot2, plot3, plot4, forMovie)

figure('visible','off')

for ct = 57:length(lagsToPlotStart)
	disp([num2str(ct)])
	lag = lagsToPlotStart(ct);
	%% If lag is > 20 days, do 10-day average %%

	lagRange = lagsToPlotEnd(ct):lagsToPlotStart(ct);
	timeAvg = bsxfun(@minus, alldays, lagRange');
	timeAvg = unique(timeAvg(:));

	timeOverlapSST = find(ismember(sst_comp.time, timeAvg));
	timeOverlapNCEP = find(ismember(hgt300_comp.time, timeAvg));

	sstComp = mean(sst_comp.anom(:,:,timeOverlapSST), 3);

	[sst_comp.LAT sst_comp.LON] = meshgrid(sst_comp.lat, sst_comp.lon);


	%%%% SST (sig) + z300 + WAF %%%%%%%
	if plot1
		%%%%% Assess significance of SST patterns %%%%%
		if length(lagRange) == 1
			sigSST = double(sstComp > sstCompSurrBnds(:,:,2,lag_range == lag) | sstComp < sstCompSurrBnds(:,:,1,lag_range == lag));
		end
		%%%%% Get composites %%%%%

		hgt300Comp = mean(hgt300_comp.anom(:,:,timeOverlapNCEP), 3);
		if lag <= 50 % on average, still in summer
			Wx300 = WAF.Wx(:, :, WAF.level == 300, WAF.leadtime == lag);
			Wy300 = WAF.Wy(:, :, WAF.level == 300, WAF.leadtime == lag);

			Wx300(abs(WAF.LAT) > 70) = NaN;
			Wy300(abs(WAF.LAT) > 70) = NaN;
		end

		%%%%% Make plot %%%%%

		clf
		if forMovie
			m_proj('miller','lat',[-10 80],'lon',[0 360]);
		else
			m_proj('miller','lat',[-10 80],'lon',[120 330]);
		end
		m_pcolorKM(sst_comp.lon, sst_comp.lat, sstComp');
		shading interp
		m_coast('patch',[.7 .7 .7],'edgecolor','none');
		m_grid('box','fancy');
		hold on
		if length(lagRange) == 1, m_plot(sst_comp.LON(sigSST > 0), sst_comp.LAT(sigSST > 0), '.k','markersize', 2); end

		dk_turq = rgb('darkturquoise');

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
		caxis([-0.4 0.4])

		clines = linspace(-80,80,17);
		m_contour(hgt300_comp.lon,hgt300_comp.lat,hgt300Comp',clines(clines>0),'linewidth',2,'col','k','linestyle','-');
		m_contour(hgt300_comp.lon,hgt300_comp.lat,hgt300Comp',clines(clines<0),'linewidth',2,'col','k','linestyle','--');
		if lag <= 50
			arrow_color = rgb('OrangeRed');
			pl = find( sqrt(Wx300.^2 + Wy300.^2) < 0.1);
			Wx300(pl) = NaN; Wy300(pl) = NaN;
			h = m_quiver(WAF.LON(1:3:end,1:3:end), WAF.LAT(1:3:end,1:3:end), Wx300(1:3:end,1:3:end), Wy300(1:3:end,1:3:end),2,'color',arrow_color, 'linewidth',2);
			if forMovie
				m_quiver(5, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
			else
				m_quiver(125, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
			end

		end
		set(gcf,'Renderer','zbuffer')
		if forMovie
			title(['Lead time: ' num2str(lag) ' days'])
			figname = ['sstSig_z300_WAF_figcount_' sprintf('%04d',ct) '.png'];
			orient landscape
			print('-dpng',[fig_folder '/' figname])
		else

			orient landscape
			print('-dpsc',[fig_folder '/sstSig_z300_WAF_lagStart' num2str(lagsToPlotStart(ct)) '_lagEnd' num2str(lagsToPlotEnd(ct)) '.ps'])
		end

	end
	%%%% SST + OA fluxes %%%%%%%
	if plot2

		timeOverlapOA = find(ismember(qnet_comp.time, timeAvg));
		qnetComp = mean(qnet_comp.anom(:,:,timeOverlapOA), 3);

		[qnet_comp.LAT qnet_comp.LON] = meshgrid(qnet_comp.lat, qnet_comp.lon);

		clf
		if forMovie
			m_proj('miller','lat',[-10 80],'lon',[0 360]);
		else
			m_proj('miller','lat',[-10 80],'lon',[120 330]);
		end
		m_pcolorKM(sst_comp.lon, sst_comp.lat, sstComp');
		shading interp
		m_coast('patch',[.7 .7 .7],'edgecolor','none');
		m_grid('box','fancy');
		hold on

		dk_turq = rgb('darkturquoise');

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
		caxis([-0.4 0.4])

		clines = -30:5:30;
		qnetComp(qnet_comp.LAT > 65) = NaN;
		% remove great lakes
		great_lakes = qnet_comp.LON >268 & qnet_comp.LON < 284 & ...
			qnet_comp.LAT > 40 & qnet_comp.LAT < 49;
		qnetComp(great_lakes) = NaN;
		m_contour(qnet_comp.lon,qnet_comp.lat,qnetComp',clines(clines>0),'linewidth',2,'col','k','linestyle','-');
		m_contour(qnet_comp.lon,qnet_comp.lat,qnetComp',clines(clines<0),'linewidth',2,'col','k','linestyle','--');

		%{
		arrow_color = rgb('OrangeRed');
		uwndComp(uwnd10m_comp.LAT > 80 | uwnd10m_comp.LAT < -10) = NaN;
		vwndComp(uwnd10m_comp.LAT > 80 | uwnd10m_comp.LAT < -10) = NaN;
		h = m_quiver(uwnd10m_comp.LON(1:3:end,1:3:end), uwnd10m_comp.LAT(1:3:end,1:3:end), ...
			uwndComp(1:3:end,1:3:end), vwndComp(1:3:end,1:3:end),2,'color',arrow_color, ...
			'linewidth',2);
		%}
		set(gcf,'Renderer','zbuffer')
		if forMovie
			title(['Lead time: ' num2str(lag) ' days'])
			figname = ['sst_qnet_figcount_' sprintf('%04d',ct) '.png'];
			orient landscape
			print('-dpng',[fig_folder '/' figname])
		else

			orient landscape
			print('-dpsc',[fig_folder '/sst_qnet_lagStart' num2str(lagsToPlotStart(ct)) '_lagEnd' num2str(lagsToPlotEnd(ct)) '.ps'])
		end
	end
	if plot3

		uwndComp = mean(uwnd10m_comp.anom(:,:,timeOverlapNCEP), 3);
		vwndComp = mean(vwnd10m_comp.anom(:,:,timeOverlapNCEP), 3);
		[uwnd10m_comp.LAT uwnd10m_comp.LON] = meshgrid(uwnd10m_comp.lat, uwnd10m_comp.lon);

		clf
		if forMovie
			m_proj('miller','lat',[-10 80],'lon',[0 360]);
		else
			m_proj('miller','lat',[-10 80],'lon',[120 330]);
		end
		m_pcolorKM(sst_comp.lon, sst_comp.lat, sstComp');
		shading interp
		m_coast('patch',[.7 .7 .7],'edgecolor','none');
		m_grid('box','fancy');
		hold on

		dk_turq = rgb('darkturquoise');

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
		caxis([-0.4 0.4])


		arrow_color = rgb('OrangeRed');
		uwndComp(uwnd10m_comp.LAT > 80 | uwnd10m_comp.LAT < -10) = NaN;
		vwndComp(uwnd10m_comp.LAT > 80 | uwnd10m_comp.LAT < -10) = NaN;
		h = m_quiver(uwnd10m_comp.LON(1:3:end,1:3:end), uwnd10m_comp.LAT(1:3:end,1:3:end), ...
			uwndComp(1:3:end,1:3:end), vwndComp(1:3:end,1:3:end),2,'color',arrow_color, ...
			'linewidth',2);
		if forMovie
			m_quiver(5, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
		else
			m_quiver(125, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
		end

		set(gcf,'Renderer','zbuffer')
		if forMovie
			title(['Lead time: ' num2str(lag) ' days'])
			figname = ['sst_surfwinds_figcount_' sprintf('%04d',ct) '.png'];
			orient landscape
			print('-dpng',[fig_folder '/' figname])
		else

			orient landscape
			print('-dpsc',[fig_folder '/sst_surfwinds_lagStart' num2str(lagsToPlotStart(ct)) '_lagEnd' num2str(lagsToPlotEnd(ct)) '.ps'])
		end
	end
	%%%% SST (sig) + z850 + WAF %%%%%%%
	if plot4
		%%%%% Assess significance of SST patterns %%%%%
		if length(lagRange) == 1
			sigSST = double(sstComp > sstCompSurrBnds(:,:,2,lag_range == lag) | sstComp < sstCompSurrBnds(:,:,1,lag_range == lag));
		end
		%%%%% Get composites %%%%%

		hgt850Comp = mean(hgt850_comp.anom(:,:,timeOverlapNCEP), 3);
		if lag <= 50 % on average, still in summer
			Wx850 = WAF.Wx(:, :, WAF.level == 850, WAF.leadtime == lag);
			Wy850 = WAF.Wy(:, :, WAF.level == 850, WAF.leadtime == lag);

			Wx850(abs(WAF.LAT) > 70) = NaN;
			Wy850(abs(WAF.LAT) > 70) = NaN;
		end

		%%%%% Make plot %%%%%

		clf
		if forMovie
			m_proj('miller','lat',[-10 80],'lon',[0 360]);
		else
			m_proj('miller','lat',[-10 80],'lon',[120 330]);
		end
		m_pcolorKM(sst_comp.lon, sst_comp.lat, sstComp');
		shading interp
		m_coast('patch',[.7 .7 .7],'edgecolor','none');
		m_grid('box','fancy');
		hold on
		if length(lagRange) == 1, m_plot(sst_comp.LON(sigSST > 0), sst_comp.LAT(sigSST > 0), '.k','markersize', 2); end

		dk_turq = rgb('darkturquoise');
 
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
		caxis([-0.4 0.4])

		clines = linspace(-40,40,17);
		m_contour(hgt850_comp.lon,hgt850_comp.lat,hgt850Comp',clines(clines>0),'linewidth',2,'col','k','linestyle','-');
		m_contour(hgt850_comp.lon,hgt850_comp.lat,hgt850Comp',clines(clines<0),'linewidth',2,'col','k','linestyle','--');
		if lag <= 50
			arrow_color = rgb('OrangeRed');
			pl = find( sqrt(Wx850.^2 + Wy850.^2) < 0.1);
			Wx850(pl) = NaN; Wy850(pl) = NaN;
			h = m_quiver(WAF.LON(1:3:end,1:3:end), WAF.LAT(1:3:end,1:3:end), Wx850(1:3:end,1:3:end), Wy850(1:3:end,1:3:end),2,'color',arrow_color, 'linewidth',2);
			if forMovie
				m_quiver(5, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
			else
				m_quiver(125, 0, 5, 0, 2,'color',arrow_color, 'linewidth',2);
			end

		end
		set(gcf,'Renderer','zbuffer')
		if forMovie
			title(['Lead time: ' num2str(lag) ' days'])
			figname = ['sstSig_z850_WAF_figcount_' sprintf('%04d',ct) '.png'];
			orient landscape
			print('-dpng',[fig_folder '/' figname])
		else

			orient landscape
			print('-dpsc',[fig_folder '/sstSig_z850_WAF_lagStart' num2str(lagsToPlotStart(ct)) '_lagEnd' num2str(lagsToPlotEnd(ct)) '.ps'])
		end

	end
end

