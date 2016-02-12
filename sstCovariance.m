function [C, C0] = sstCovariance(GHCND_clustered, sst_comp, alldays, lag_range, domainlatlims,...
	domainlonlims, yrs, nyrs, summer_doy, leaveNout, excludeYear, makeplot, fig_folder)
% Calculate the spatial covariance between SST anomalies on a given day, and the identified
% Pacific SST pattern based on the subsequent observation of hot days.
% Analysis is done via 'leave one out' such that, when calculating the covariance for a
% given year, the data from that year is not used to produce the composite.
   
[LAT LON] = meshgrid(sst_comp.lat, sst_comp.lon);

% get data in the desired region
domainInd = LAT > min(domainlatlims) & LAT < max(domainlatlims) & ...
	LON > min(domainlonlims) & LON < max(domainlonlims);
sstVec = reshape(sst_comp.anom,[length(sst_comp.lon)*length(sst_comp.lat) length(sst_comp.time)]);
sstVec = sstVec(domainInd,:);
 
lonPEP = sst_comp.lon(sst_comp.lon > min(domainlonlims) & sst_comp.lon < max(domainlonlims));
latPEP = sst_comp.lat(sst_comp.lat > min(domainlatlims) & sst_comp.lat < max(domainlatlims));


halfWindow = floor(leaveNout/2);

% pull out land, if present
pl = ~isnan(sstVec(:, 1));
sstVec = sstVec(pl, :);
 
% summer time
summerDates = GHCND_clustered.time(ismember(GHCND_clustered.doy, summer_doy));

% dates of hot days (all)
hotDays = GHCND_clustered.time(ismember(GHCND_clustered.time, alldays));

C = NaN(length(lag_range), length(summerDates));

for jj = 1:length(lag_range)
	lag = lag_range(jj);
	PEPsave = NaN(sum(domainInd(:)), length(yrs));
	for ii = 1:nyrs
		yrExclude = [(yrs(ii)-halfWindow):(yrs(ii)+halfWindow), excludeYear]; 
		
		hotDaysUse = hotDays(~ismember(year(hotDays), yrExclude));
		% composite on summers except the one being predicted
		idxSST = find(ismember(sst_comp.time, hotDaysUse)); % sst on hot days except for the summer being predicted
		idxLagged = idxSST - lag; % ssts lead hot days
		idxLagged(idxLagged < 1) = [];
		idxLagged(idxLagged > size(sstVec,2)) = [];
		% composite and remove spatial mean for covariance calculation
		sstComp = detrend(mean(sstVec(:,idxLagged),2),'constant');
		
		if makeplot(jj) 
			PEPsave(:,ii) = sstComp;
		end

		% get spatial covariance for summer being predicted
		% Pull out locations that lead the summer days
		idxPredict = find(ismember(year(sst_comp.time), yrs(ii)) & ...
			ismember(sst_comp.doy, summer_doy)) - lag;
		pl1 = idxPredict < 1; 
		pl2 = idxPredict > size(sstVec, 2);
		idxPredict(pl1) = [];
		idxPredict(pl2) = [];
		sstAnom = detrend(sstVec(:,idxPredict),'constant');
		numvals = sum(~isnan(sstAnom), 1);

		idxUse = find(ismember(year(summerDates), yrs(ii)));
		idxUse(pl1 | pl2) = [];
		C(jj, idxUse) = sstComp'*sstAnom./(numvals - 1);

	end
	if makeplot(jj)
		disp(['Making plot for lead time of ' num2str(lag) ' days'])
		Pstd = std(PEPsave,[],2);
		Pstd = reshape(Pstd,[numel(lonPEP) numel(latPEP)]);

		clf
		m_proj('miller','lat',[-10 80],'lon',[120 330]);
		hold on
		m_pcolorKM(lonPEP, latPEP, Pstd');
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
		cmap = flipud(lbmap(17, 'brownblue'));
		cmap(9,:) = [];
		colormap(cmap)
		h = colorbar;
		caxis([0 0.1])
		
		set(gcf,'Renderer','zbuffer')
		figname = ['PEPstd_lag' num2str(lag) '.ps'];
		orient landscape
		print('-dpsc',[fig_folder '/' figname]);
	end
end

C0 = C;
% normalize C by the standard deviation for each lag
C = bsxfun(@times, C, 1./std(C,[],2));

