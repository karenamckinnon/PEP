function precip_comp = precip_index(precip_comp, GHCND_clustered, avgtime0, getspi, getSPI)
% Calculate average precipitation anomaly in cluster, as well as cluster-wide
% and gridbox SPI
% avgtime = number of days to use for calculating SPI

[LAT LON] = meshgrid(precip_comp.lat, precip_comp.lon);
nlat = length(precip_comp.lat);
nlon = length(precip_comp.lon);
ntime = length(precip_comp.time);

%% CALCULATE AVERAGE PRECIPITATION ANOMALY IN THE CLUSTER %%

% get polygon that roughly encloses east coast points
poly_indices = boundary(GHCND_clustered.location(:,1), GHCND_clustered.location(:,2));
boundary_pts = [GHCND_clustered.location(poly_indices,1) GHCND_clustered.location(poly_indices,2)];
% switch to [0 360]
boundary_pts(:,1) = 360 + boundary_pts(:,1);
points_in_eastern_us = inpolygon(LON(:), LAT(:), boundary_pts(:,1), boundary_pts(:,2));

% calculate average precipitation anomalies in the region
Zmat = reshape(double(precip_comp.anom),[nlon*nlat ntime]);
avgPrecip = nanmean(Zmat(points_in_eastern_us,:),1);
precip_comp.avgPrecip = avgPrecip;

% and average precip, generally
Zmat = reshape(double(precip_comp.data0),[nlon*nlat ntime]);
avgPrecipTotal = nanmean(Zmat(points_in_eastern_us,:),1);
precip_comp.avgPrecipTotal = avgPrecipTotal;

precip_comp.doy = precip_comp.time - datenum(year(precip_comp.time), 1, 1) + 1;

%% CALCULATE THE SPI FOR THE REGION %%

% uses precipitation data (not anomalies)
Zmat = reshape(double(precip_comp.data0),[nlon*nlat ntime]);
avgPrecip = nanmean(Zmat(points_in_eastern_us,:),1);
yrs = unique(year(precip_comp.time));

for counter = 1:length(avgtime0)
	disp(['Calculating SPI for averaging time of ' num2str(avgtime0(counter)) ' days'])
	tic;
	avgtime = avgtime0(counter);
	if getspi(counter)
		ts = avgPrecip; % average for region
		spi = NaN(size(ts));
		unique_doy = unique(precip_comp.doy);
		for doy = 1:max(unique_doy)

			doy_range = (doy-avgtime + 1):(doy); % calculate over 'avgtime' preceding
			ndoy = length(doy_range);
			doy_range(doy_range < 1) = doy_range(doy_range < 1) + 365;
			if doy ~= 366
				doy_range(doy_range > 365) = doy_range(doy_range > 365) - 365;
				daily_vals = ts(ismember(precip_comp.doy, doy_range));
			else
				doy_range(doy_range > 366) = doy_range(doy_range > 366) - 366;
				idx1 = strfind(precip_comp.doy',doy_range);
				idx2 = find(precip_comp.doy == 366);
				daily_vals = [];
				for kk = 1:length(idx1)
					daily_vals = [daily_vals ts(idx1(kk):idx2(kk))];
				end
			end

			index = find(ismember(precip_comp.doy, doy));

			% catch if there are not the same # of days per year 
			% can happen on edges of dataset
			if mod(length(daily_vals)/ndoy,1) ~= 0 
				doy_use = precip_comp.doy(ismember(precip_comp.doy, doy_range));
				indices = [find(doy_use == doy_range(1),1,'first') find(doy_use == doy_range(end),1,'last')];
				daily_vals = daily_vals(indices(1):indices(2));
				if doy_use(1) ~= doy_range(1),
					index(1) = [];
				end
			end
			 
			% get the sum of the daily values
			summedVals = mean(reshape(daily_vals, [ndoy length(daily_vals)/ndoy]), 1);
			
			% fit gamma
			gamma_params = gamfit(summedVals);

			gamCdfVal = gamcdf(summedVals, gamma_params(1), gamma_params(2));

			% SPI is the map from the gamma to the normal
			spi(index) = norminv(gamCdfVal);
		end

		spiname = ['spi' num2str(avgtime) '']
		precip_comp.(spiname) = spi;
	end

	%% CALCULATE THE SPI FOR EACH GRIDBOX (TAKES A LONG TIME) %%
	if getSPI(counter)
		rng('default');
		spi_input.time = precip_comp.time;
		spi_input.avgtime = avgtime;
		spi_hash = DataHash(spi_input);
		savename =['/n/huybers_lab/common/ghcnd/analysis/composites/cache/SPI.' spi_hash '.nc'];
		SPIname = ['SPI' num2str(avgtime) '']

		if exist(savename,'file')
			disp(['Loading SPI for an averaging time of ' num2str(avgtime) ''])
			SPI = ncread(savename,'SPI');
			precip_comp.(SPIname) = SPI;
		else
			disp(['Need to calculate SPI for an averaging time of ' num2str(avgtime) ''])
			loc = find(~isnan(Zmat(:,1)));
			SPI = NaN(size(Zmat));
			for ct = 1:length(loc)
				% disp(['Calculating SPI for location ' num2str(ct) ' of ' num2str(length(loc)) ''])
				ts = Zmat(loc(ct),:); % time series at a given gridbox
				spi = NaN(size(ts));
				for doy = 1:max(unique_doy)
					doy_range = (doy-avgtime):(doy); % calculate over the month preceding
					ndoy = length(doy_range);
					doy_range(doy_range < 1) = doy_range(doy_range < 1) + 365;
					if doy ~= 366
						doy_range(doy_range > 365) = doy_range(doy_range > 365) - 365;
						daily_vals = ts(ismember(precip_comp.doy, doy_range));
					else
						doy_range(doy_range > 366) = doy_range(doy_range > 366) - 366;
						idx1 = strfind(precip_comp.doy',doy_range);
						idx2 = find(precip_comp.doy == 366);
						daily_vals = [];
						for kk = 1:length(idx1)
							daily_vals = [daily_vals ts(idx1(kk):idx2(kk))];
						end
					end

					index = find(ismember(precip_comp.doy, doy));

					if mod(length(daily_vals)/ndoy,1) ~= 0
						doy_use = precip_comp.doy(ismember(precip_comp.doy, doy_range));
						indices = [find(doy_use == doy_range(1),1,'first') find(doy_use == doy_range(end),1,'last')];
						daily_vals = daily_vals(indices(1):indices(2));
						if doy_use(1) ~= doy_range(1),
							index(1) = [];
						end
					end
					summedVals = mean(reshape(daily_vals, [ndoy length(daily_vals)/ndoy]), 1);
					if any(summedVals == 0), % add a tiny amount of rain so parameters can be estimated
						loc0 = find(summedVals == 0);
						summedVals(loc0) = abs(1e-10*randn(size(loc0)));
					end
					% fit gamma
					gamma_params = gamfit(summedVals);

					gamCdfVal = gamcdf(summedVals, gamma_params(1), gamma_params(2));

					spi(index) = norminv(gamCdfVal);
				end
				SPI(loc(ct),:) = spi;
			end
			SPI = reshape(SPI,[nlon nlat ntime]);
			dimnames = {'lon','lat','time'};
			dims = {precip_comp.lon, precip_comp.lat, precip_comp.time};
			varnames = {'SPI'};
			vars = {SPI};
			saveVars2NetCDF(savename, dimnames, dims, varnames, vars)
			precip_comp.(SPIname) = SPI;
		end
	end
	toc

end



