% main.m
% Wrapper function for all analyses in McKinnon, Rhines, Tingley, and Huybers
% "Long-lead predictions of Eastern US hot days from Pacific sea surface temperatures"
% Necessary data files
% (1) .mat files for all GHCND stations
% (2) .nc files for other SST, z300, z850, atmos-ocean fluxes, surface winds

%%%%%%%% User inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put your own data directory here
dataDir = ''; 
% where to save figures?
fig_folder = ''; 
% where are the GHCND matfiles
ghcnd_matpath = '';
% where to save cached files?
cacheDir = '';

% all of the below fields are anomalies from the climatology
% anomalies are calculated using three annual harmonics using scripts in ncl
% one example script is provided for 10m uwnd (uwnd10m_anomalies.ncl) and all other calculations are analagous
sst_fname = [dataDir '/OISST_V2/SST/sst.anomalies.1981.2015.nc'];
hgt300_fname = [dataDir '/NCEP_R2/hgt/hgt300.anomalies.1979.2015.nc'];
hgt850_fname = [dataDir '/NCEP_R2/hgt/hgt850.anomalies.1979.2015.nc'];
qnet_fname = [dataDir '/OAFlux/qnet.anomalies.1985.2009.nc'];
surfU_fname = [dataDir '/NCEP_R2/uwnd/10m/uwnd.10m.anomalies.1979.2015.nc'];
surfV_fname = [dataDir '/NCEP_R2/vwnd/10m/vwnd.10m.anomalies.1979.2015.nc'];
hadsst_fname = [dataDir '/HadSST3/HadSST.3.1.1.0.median.nc'];

% data acquired via get_cpc_data.m
% seasonal cycle not removed at this point
precip_fname = [dataDir 'cpc-precip/precip.V1.0.1979.2015.nc'];

%%%%%%%% Parameters for the analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_begin = 1982; 
yr_end = 2015;
varname = 'TMAX'; 
baseline_percentile = 95; % threshold for hot weather
above_percentile = 1; % hot weather is above the 95th percentile 
summer_doy = 175:234; % hottest 60 days of summer
station_latlims = [24 50]; % domain for weather stations
station_lonlims = [-125 -60];
domainlatlims = [20 50]; % domain for PEP
domainlonlims = [145 230];
years_to_average = 0; % 0: = remove linear trend; val > 0 = remove running mean
lag_range = -10:60; % how many lags to calculate predictive skill for
predictandName = 'alldays'; % 'alldays' to predict hot days, 'startDate' to predict starting date of heat events
leaveNout = 1; % cross-validation
cluster_to_use = 1; % use biggest cluster
maxclust = 10; % out of 10 total clusters
cutoffPerc = 0:10:100; % thresholds for predictions (in percentile space)

figure('visible','off')

%%%%%%%% Get GHCND data, and do QC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get GHCND data
GHCND = get_ghcnd(yr_begin,yr_end, varname, ghcnd_matpath, station_latlims, station_lonlims);

%% Do QC: Require 80% of days to be present for 80% of the summer seasons (JJA)
GHCND_QC = do_qc(GHCND);

% get station areas for all stations
[station_areas, stn_polylon, stn_polylat] = tessellate(GHCND_QC, cacheDir);

%% Remove Canadian stations
canada = ismember(GHCND_QC.id(:,1),'C');
GHCND_QC.data(canada,:) = [];
GHCND_QC.location(canada,:) = [];
GHCND_QC.id(canada,:) = [];
station_areas(canada) = [];

GHCND_QC.doy = GHCND_QC.time - datenum(year(GHCND_QC.time), 1, 1) + 1;

summerDays = ismember(GHCND_QC.doy, summer_doy);
summerLength = length(summer_doy);
ntime = length(GHCND_QC.time);
nstations = size(GHCND_QC.id,1);
yrs = unique(year(GHCND_QC.time));
nyrs = length(yrs);

%% Remove climatology at each station using first three harmonics
GHCND_QC = remove_climo_ghcnd(GHCND_QC, years_to_average, cacheDir);

%%%%%%%% Cluster stations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cluster stations based on whether they are above the 95th percentile
cutoffHeat = prctile(GHCND_QC.anomTotal(:, summerDays), 95, 2);
GHCND_QC.anom = bsxfun(@minus, GHCND_QC.anomTotal, cutoffHeat);
plotDendro = 0; % plot dendrogram? Takes a long time
GHCND_clustered0 = cluster_on_extremes(GHCND_QC, maxclust, summer_doy, plotDendro, ...
	station_latlims, station_lonlims, fig_folder, cacheDir);

clusterInd = GHCND_clustered0.clusterid == cluster_to_use;
eastInd = GHCND_clustered0.location(:,1) > -105; % lone station farther west
GHCND_clustered = subset(GHCND_clustered0, clusterInd & eastInd);
station_areas0 = station_areas;
station_areas = station_areas0(clusterInd & eastInd);
areaThreshold = prctile(station_areas, 98); % cap stations area
station_areas(station_areas > areaThreshold) = areaThreshold;

%%%%%%%% Calcuate T95 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate T95 as the spatial 95th percentile in temperature anomalies
% T95prime has zero mean
[T95prime, T95] = findHotRegion(GHCND_clustered, baseline_percentile, ...
	station_areas, summerDays);

%% Define heat periods as periods where T95prime > 0 for at least
% two consecutive days, and discard periods where there is less than
% 'cutofflength' days between the end of one and the beginning of the next
cutofflength = 3;
[startDate, endDate, totalTemp, alldays] = ...
	findHotDays(GHCND_clustered, cutofflength, T95prime, summerDays, ...
	yrs, summer_doy);
	
%%%%%%%% Calcuate predictors (SST covariance) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the covariance between SST anoms and PEP to use as predictor
% Called 'C' in this script, and 'P' in the text of the paper
rng('default')
clear cache_input
cache_input.sst_fname = sst_fname;
cache_input.alldays = alldays;
cache_input.lag_range = lag_range;
cache_input.domainlatlims = domainlatlims;
cache_input.domainlonlims = domainlonlims;
cache_input.yrs = yrs;
cache_input.leaveNout = leaveNout;
cache_input.excludeYear = []; % remove one year from the analysis? 
covHash = DataHash(cache_input);
matname = [cacheDir '/cache/cov' covHash '.mat'];
if ~exist(matname, 'file')
	disp('Calculating PEP index...')
	%% Get met fields, and remove the running mean
	%% Get SST data
	sst_comp = fetch_composite(sst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
	sst_comp = removeTrend(sst_comp);
	sst_comp.doy = sst_comp.time - datenum(year(sst_comp.time), 1, 1) + 1;

	makeplot = ismember(lag_range,[50 40 30 20 15 10 5 0]);
	makeplot = zeros(size(makeplot));
	[C, C0] = sstCovariance(GHCND_clustered, sst_comp, alldays, lag_range, domainlatlims,...
		domainlonlims, yrs, nyrs, summer_doy, leaveNout, cache_input.excludeYear, makeplot, ...
		fig_folder);
	timestamp = datestr(now);
	save(matname, 'C', 'C0','timestamp');
else
	load(matname)
end

%% ATLANTIC
domainlatlimsATL = [37.5 52.5]; % from Donat et al (2015)
domainlonlimsATL = [302.5 317.5];

rng('default')
clear cache_input
cache_input.sst_fname = sst_fname;
cache_input.alldays = alldays;
cache_input.lag_range = lag_range;
cache_input.yrs = yrs;
cache_input.domainlatlims = domainlatlimsATL;
cache_input.domainlonlims = domainlonlimsATL;
cache_input.leaveNout = leaveNout;
cache_input.excludeYear = [];
covHash = DataHash(cache_input);
matname = [cacheDir '/cov' covHash '.mat'];
if ~exist(matname, 'file')
	if ~exist('sst_comp','var')
		sst_comp = fetch_composite(sst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
		sst_comp = removeTrend(sst_comp);
		sst_comp.doy = sst_comp.time - datenum(year(sst_comp.time), 1, 1) + 1;
	end
	[CATL, C0ATL] = sstCovariance(GHCND_clustered, sst_comp, alldays, lag_range, domainlatlimsATL,...
		domainlonlimsATL, yrs, nyrs, summer_doy, leaveNout, cache_input.excludeYear);
	timestamp = datestr(now);
	save(matname, 'CATL', 'C0ATL','timestamp');
else
	load(matname)
end

%% TROPICS, NO 1988

domainlatlimsNINO = [-5 5]; % Nino3.4 bounds
domainlonlimsNINO = [190 240];

rng('default')
clear cache_input
cache_input.sst_fname = sst_fname;
cache_input.alldays = alldays;
cache_input.lag_range = lag_range;
cache_input.yrs = yrs;
cache_input.domainlatlims = domainlatlimsNINO;
cache_input.domainlonlims = domainlonlimsNINO;
cache_input.leaveNout = leaveNout;
cache_input.excludeYear = 1988;
covHash = DataHash(cache_input);
matname = [cacheDir '/cov' covHash '.mat'];
if ~exist(matname, 'file')
	%% Get SST data
	if ~exist('sst_comp','var')
		sst_comp = fetch_composite(sst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
		sst_comp = removeTrend(sst_comp);
		sst_comp.doy = sst_comp.time - datenum(year(sst_comp.time), 1, 1) + 1;
	end
	[CNINO_NO88, C0NINO_NO88] = sstCovariance(GHCND_clustered, sst_comp, alldays, lag_range, domainlatlimsNINO,...
		domainlonlimsNINO, yrs, nyrs, summer_doy, leaveNout, cache_input.excludeYear);
	timestamp = datestr(now);
	save(matname, 'CNINO_NO88', 'C0NINO_NO88','timestamp');
else
	load(matname)
end

%% TROPICS, WITH 1988
rng('default')
cache_input.excludeYear = [];
covHash = DataHash(cache_input);
matname = [cacheDir '/cov' covHash '.mat'];
if ~exist(matname, 'file')
	%% Get SST data
	if ~exist('sst_comp','var')
		sst_comp = fetch_composite(sst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
		sst_comp = removeTrend(sst_comp);
		sst_comp.doy = sst_comp.time - datenum(year(sst_comp.time), 1, 1) + 1;
	end
	[CNINO_NO88, C0NINO_NO88] = sstCovariance(GHCND_clustered, sst_comp, alldays, lag_range, domainlatlimsNINO,...
		domainlonlimsNINO, yrs, nyrs, summer_doy, leaveNout, cache_input.excludeYear);
	timestamp = datestr(now);
	save(matname, 'CNINO', 'C0NINO','timestamp');
else
	load(matname)
end

%%%%%%%% Assess skill with ROC scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate null bounds for the ROC score using PEP as a predictor by bootstrapping
rng('default')
clear cache_input
summerTime = GHCND_clustered.time(ismember(GHCND_clustered.doy, summer_doy));
nsurr = 1e4;
rocScoreSurr = NaN(nsurr, length(lag_range));
makeplots = 0;
cache_input.nsurr = nsurr;
cache_input.alldays = alldays;
cache_input.lag_range = lag_range;
cache_input.metric = predictandName;
cache_input.C = C;
ROCsurrHash = DataHash(cache_input);
matname = [cacheDir '/ROCsurr' ROCsurrHash '.mat'];
if exist(matname,'file')
	load(matname)
else
	for ct = 1:nsurr
		disp([num2str(ct)])
		idxSurr = datasample(1:34,34);
		if strcmp(cache_input.metric,'alldays')
			hotMat = ismember(GHCND_clustered.time, alldays);
		elseif strcmp(cache_input.metric,'startDate')
			hotMat = ismember(GHCND_clustered.time, startDate);
		end 
		hotMat = reshape(hotMat(ismember(GHCND_clustered.doy, summer_doy)),[60 34]);
		hotMatSurr = hotMat(:, idxSurr);
		predictandSurr = summerTime(hotMatSurr(:));
		[rocScoreSurr(ct, :)] = ...
			calcROC(GHCND_clustered, C, [], lag_range, predictandSurr, ...
			summer_doy, makeplots, fig_folder, 'pacific', cutoffPerc);

	end
	save(matname,'rocScoreSurr');
end

%% Make plots of ROC scores
lagsToPlot = [20 30 40 50];
predictand = eval(predictandName);
[roc_score,TPR, FPR, cutoffs] = pepROC(lag_range, lagsToPlot, C, 'PEP', summer_doy, predictand, predictandName, ...
	GHCND_clustered, fig_folder);

%% significance of ROC scores as compared to surrogate
for ct = 1:length(lag_range)
	pvalPEP(ct) = length(find(rocScoreSurr(:,ct) > roc_score(ct)))/nsurr;
end 

%% ROC scores for Atlantic and tropics as predictor
[roc_scoreATL] = pepROC(lag_range, lagsToPlot, CATL, 'ATL', summer_doy, predictand, predictandName, ...
	GHCND_clustered, fig_folder);

% For predictions from tropics, compare the model with and without 1988
% need to do this because 1988 strongly skews predictions (i.e. makes them especially bad)
% because it was both a strong La Nina and a very hot year in the Eastern US
% however, La Ninas are not *generally* associated with hot Eastern US summers
Cmat = NaN(2, size(CNINO, 1), size(CNINO, 2)); % put both with and without 1988 in same array
Cmat(1, :, :) = CNINO;
Cmat(2, :, :) = CNINO_NO88;
[roc_scoreNINO] = pepROC(lag_range, lagsToPlot, Cmat, 'NINO', summer_doy, predictand, predictandName, ...
	GHCND_clustered, fig_folder);

%% Calculate and plot ROC scores per station
roc_scoreStation = plotROCstations(GHCND_clustered, cutofflength, C, lag_range, lagsToPlot, ...
	fig_folder, summerDays, yrs, summer_doy, cone_ratio, [predictandName 'Station']);

%%%%%%%% Do similar analysis using SPI as predictor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Precipitation analysis using CPC data, and mapping to SPI
rng('default')
clear cache_input
avgtime0 = [15 30 60 90 180 270];
getspi = ones(size(avgtime0)); % get SPI for region as a whole
getSPI = [0 1 1 0 0 0]; % get SPI for individual gridboxes. Already calculated, but takes a long time to load
timePrecip = datenum(1981,9,1):datenum(2015,8,31);
cache_input.precip_fname = precip_fname;
cache_input.avgtime0 = avgtime0;
cache_input.getSPI = getSPI;
cache_input.getspi = getspi;
cache_input.timePrecip = timePrecip;
precipHash = DataHash(cache_input);
matname = [cacheDir '/precip' precipHash '.mat'];
if ~exist(matname,'file')
	precip_comp = fetch_composite(precip_fname, timePrecip, 'precip');
	precip_comp = process_precip(precip_comp, precip_fname); % some issues with data file
	precip_comp = remove_climatology(precip_comp, 'precip');
	precip_comp.data0 = precip_comp.data;
	precip_comp.data = precip_comp.anom; % change the names here to read into removeTrend
	precip_comp = removeTrend(precip_comp);
	disp('Got precip data')
	precip_comp.doy = precip_comp.time - datenum(year(precip_comp.time), 1, 1) + 1;

	precip_comp = precip_index(precip_comp, GHCND_clustered, avgtime0, getspi, getSPI);
	precip_comp0 = precip_comp;
	precip_comp = rmfield(precip_comp,{'data','time_nc','anom','data0'});
	save(matname,'-v7.3','precip_comp'); % big file! 
else 
	load(matname)
end

%% Make ROC plots for precip
varsToPlot = {precip_comp.spi30, precip_comp.spi60, precip_comp.spi90, precip_comp.spi180};
names = {'SPI30','SPI60','SPI90','SPI180'};
summerPrecip = find(ismember(precip_comp.doy, summer_doy));
for jj = 1:length(varsToPlot)
	fld = varsToPlot{jj};
	CPrecip = NaN(length(lag_range), length(summerPrecip));
	for kk = 1:length(lag_range)
		try
			CPrecip(kk,:) = fld(summerPrecip - lag_range(kk));
		end
	end
	CPrecip = -bsxfun(@times, CPrecip, 1./nanstd(CPrecip,[],2));
	makeplots = 0;
	
	[roc_scorePRECIP] = pepROC(lag_range, lagsToPlot, CPrecip, names{jj}, summer_doy, predictand, predictandName, ...
		GHCND_clustered, fig_folder);
end

% get values for ROC, TPR, FPR, etc. for all of the SPIs
[roc_scorePrecip, TPRprecip, FPRprecip, cutoffs, CPrecip] = precipROC(lag_range, [], ...
	varsToPlot, names, precip_comp, summer_doy, predictand, ...
	predictandName, GHCND_clustered, fig_folder);
	
%% Calculate null bounds for precip
rng('default')
clear cache_input
summerTime = GHCND_clustered.time(ismember(GHCND_clustered.doy, summer_doy));
nsurr = 1e4;
rocScoreSurrPrecip = NaN(nsurr, length(lag_range));
makeplots = 0;
cache_input.nsurr = nsurr;
cache_input.alldays = alldays;
cache_input.lag_range = lag_range;
cache_input.metric = predictandName;
cache_input.C = CPrecip;
ROCsurrHash = DataHash(cache_input);
matname = [cacheDir '/ROCsurr' ROCsurrHash '.mat'];
if exist(matname,'file')
	load(matname)
else
	for ct = 1:nsurr
		disp([num2str(ct)])
		idxSurr = datasample(1:34,34);
		if strcmp(cache_input.metric,'alldays')
			hotMat = ismember(GHCND_clustered.time, alldays);
		elseif strcmp(cache_input.metric,'startDate')
			hotMat = ismember(GHCND_clustered.time, startDate);
		end
		hotMat = reshape(hotMat(ismember(GHCND_clustered.doy, summer_doy)),[60 34]);
		hotMatSurr = hotMat(:, idxSurr);
		predictandSurr = summerTime(hotMatSurr(:));
		[rocScoreSurrPrecip(ct, :)] = ...
			calcROC(GHCND_clustered, cache_input.C, [], lag_range, predictandSurr, ...
			summer_doy, makeplots, fig_folder, 'pacific', cutoffPerc);

	end
	save(matname,'rocScoreSurrPrecip');
end

%% Get significance for precip predictors
for ct = 1:length(lag_range)
	pvalPrecip(ct) = length(find(rocScoreSurrPrecip(:,ct) > roc_scorePrecip(ct)))/nsurr;
end

%%%%%%%% Make figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load files for figures
if ~exist('sst_comp','var')
	sst_comp = fetch_composite(sst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
	sst_comp = removeTrend(sst_comp);
	sst_comp.doy = sst_comp.time - datenum(year(sst_comp.time), 1, 1) + 1;
end

%% Get z300 data
hgt300_comp = fetch_composite(hgt300_fname, sst_comp.time, 'hgt');
hgt300_comp = removeTrend(hgt300_comp);

%% Get z850 data
hgt850_comp = fetch_composite(hgt850_fname, sst_comp.time, 'hgt');
hgt850_comp = removeTrend(hgt850_comp);

%% Get OA fluxes
qnet_comp = fetch_composite(qnet_fname, sst_comp.time, 'qnet');
qnet_comp = removeTrend(qnet_comp);

%% Get surface winds
uwnd10m_comp = fetch_composite(surfU_fname, sst_comp.time, 'uwnd');
uwnd10m_comp = removeTrend(uwnd10m_comp);

vwnd10m_comp = fetch_composite(surfV_fname, sst_comp.time, 'vwnd');
vwnd10m_comp = removeTrend(vwnd10m_comp);

%% Get WAF for all hot days
% this calls NCL
% mean state: JJA
rng('default')
clear cache_input
nc_units = 'hours since 1800-1-1 00:00:0.0';
alldaysNC = 24*(alldays - datenum(1800,1,1));
cache_input.lag_range = -10:100;
cache_input.days = alldaysNC;
WAFhash = DataHash(cache_input);
savename = ['hotdays.' WAFhash '.csv'];
ncsavename = [cacheDir '/takaya.flux.composites.3d.' WAFhash '.nc'];
if ~exist(ncsavename, 'file')
	dlmwrite([cacheDir '/' savename],alldaysNC,'precision','%7i');
	system(['ncl lagMin=' num2str(min(lag_range)) ' lagMax=' num2str(max(lag_range)) ' ''ncName="' ncsavename '"'' ''hotDaysFname="' savename '"'' compute_W_comp3d.ncl']);
end

WAF.lat = ncread(ncsavename,'lat');
WAF.lon = ncread(ncsavename,'lon');
WAF.level = ncread(ncsavename,'level');
WAF.Wx = ncread(ncsavename,'Wx');
WAF.Wy = ncread(ncsavename,'Wy');
WAF.Wz = ncread(ncsavename,'Wz');
WAF.leadtime = ncread(ncsavename,'time');
[WAF.LAT WAF.LON] = meshgrid(WAF.lat, WAF.lon);

%% Create surrogate SST pattern based on bootstrapping years to assess significance of anomalies
nsurr = 1e3;
prctilesSig = [5 95];
[sstCompSurrBnds] = sstSurrogates(nsurr, sst_comp, summer_doy, alldays, ...
	summerLength, nyrs, prctilesSig, lag_range);
	
% Plot composite pattern of z300, SST, precip (Fig. S14)
precip_cutoff = 0.5; % define wet/dry as >/< 0.5 SPI
precipCz300(C, lag_range, GHCND_clustered, precip_comp, hgt300_comp, ...
	sst_comp, summerDays, precip_cutoff, fig_folder)

%% Plot composite maps for different lead times
forMovie = 1; % movie plots .png; no movie plots .ps
plot1 = 1; % SST + z300 + WAF
plot2 = 1; % SST + OAflux
plot3 = 1; % SST + surfwinds
plot4 = 1; % SST + z850 + WAF
lagsToPlotStart = [50:-1:-10];
lagsToPlotEnd =   [50:-1:-10];

makeMaps(lag_range, lagsToPlotStart, lagsToPlotEnd, sst_comp, hgt300_comp, ...
	hgt850_comp, uwnd10m_comp, vwnd10m_comp, ...
	qnet_comp, alldays, WAF, fig_folder, ...
	domainlatlims, domainlonlims, sstCompSurrBnds, plot1, plot2, plot3, plot4, forMovie)

%% plot the intra seasonal variability for each year
predictand = eval(predictandName);
yearByYearFigs(C, T95, GHCND_clustered, summerTime, yrs, summer_doy, ...
	lag_range, fig_folder, predictand)

% plot case studies
date_to_plot = [startDate(endDate - startDate > 7)];
lags_to_plot = [50 40 30 20 15 10 5 0];
lags_to_plot = [20];
plot_case_studies(C, GHCND_clustered, ...
	sst_comp, hgt300_comp, lags_to_plot, lag_range, summer_doy, ...
	domainlatlims, domainlonlims, date_to_plot(2), fig_folder)

%%%%%%%% HadSST3 analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% train on same period as before
hadsst_comp = fetch_composite(hadsst_fname, datenum(1981,9,1):datenum(2015,8,31), 'sst');
hadsst_comp = removeClimoHadSST(hadsst_comp, 'sst');
dataMat = reshape(hadsst_comp.anom,[size(hadsst_comp.anom,1)*size(hadsst_comp.anom,2) size(hadsst_comp.anom,3)]);
hadsst_comp.anom = reshape(nandetrend(dataMat')',size(hadsst_comp.anom));

% use to predict earlier time period
years_test = [1950:1981];

%% Get Hadley data for longer time period
GHCND_had = get_ghcnd(years_test(1),years_test(end), varname, ghcnd_matpath, station_latlims, station_lonlims);
%% Do QC: Require 80% of days to be present for 80% of the summer seasons (JJA)
GHCND_had_QC = do_qc(GHCND_had);

% get station areas for all stations
[station_areas2] = tessellate(GHCND_had_QC);

% get stations that are also present in the main analysis
%% Pick out subset of stations that passes QC and is in the originally defined cluster
[GHCND_clustered_had, overlapIdx] = match_stations(GHCND_had_QC, GHCND_clustered);
station_areas2 = station_areas2(overlapIdx);

GHCND_clustered_had.doy = GHCND_clustered_had.time - datenum(year(GHCND_clustered_had.time), 1, 1) + 1;

%% Remove climatology at each station using first three harmonics
GHCND_clustered_had = remove_climo_ghcnd(GHCND_clustered_had, years_to_average);

areaThreshold = prctile(station_areas2, 98);
station_areas2(station_areas2 > areaThreshold) = areaThreshold;

summerDaysHad = ismember(GHCND_clustered_had.doy, summer_doy);
ntimeHad = length(GHCND_clustered_had.time);
nstationsHad = size(GHCND_clustered_had.id,1);
yrsHad = unique(year(GHCND_clustered_had.time));
nyrsHad = length(yrsHad);

% Get T95
[T95primeHad] = findHotRegion(GHCND_clustered_had, baseline_percentile, ...
	station_areas2, summerDaysHad);

%% Define heat periods as periods where T95prime > 0 for at least
% two consecutive days, and discard periods where there is less than
% 'cutofflength' days between the end of one and the beginning of the next

cutofflength = 3;

[startDateHad, endDateHad, totalTempHad, alldaysHad] = ...
	findHotDays(GHCND_clustered_had, cutofflength, T95primeHad, summerDaysHad, ...
	yrsHad, summer_doy);

hadsst_comp2 = fetch_composite(hadsst_fname, datenum(years_test(1)-1,1,1):datenum(years_test(end),12,31), 'sst');
hadsst_comp2 = removeClimoHadSST(hadsst_comp2, 'sst');
dataMat = reshape(hadsst_comp2.anom,[size(hadsst_comp2.anom,1)*size(hadsst_comp2.anom,2) size(hadsst_comp2.anom,3)]);
hadsst_comp2.anom = reshape(nandetrend(dataMat')',size(hadsst_comp2.anom));

% switch longitude to 0,360
lon = hadsst_comp.lon;
lon(lon < 0) = lon(lon < 0) + 360;
[lon, idx] = sort(lon, 'ascend');
hadsst_comp2.anom = hadsst_comp2.anom(idx,:,:);
hadsst_comp2.lon = hadsst_comp2.lon(idx);
hadsst_comp.anom = hadsst_comp.anom(idx,:,:);
hadsst_comp.lon = hadsst_comp.lon(idx);

nlon = numel(hadsst_comp.lon);
nlat = numel(hadsst_comp.lat);
[LAT LON] = meshgrid(hadsst_comp.lat,hadsst_comp.lon);
domain_ind = find(LAT > min(domainlatlims) & LAT < max(domainlatlims) &...
    LON > min(domainlonlims) & LON < max(domainlonlims));

Zmat = reshape(hadsst_comp.anom, [nlon*nlat size(hadsst_comp.anom, 3)]);
Zmat = Zmat(domain_ind, :); % pull out domain of interest
Zmat2 = reshape(hadsst_comp2.anom, [nlon*nlat size(hadsst_comp2.anom, 3)]);
Zmat2 = Zmat2(domain_ind, :); % pull out domain of interest

lagRangeHad = [1 2]; % months
% do predictions but only at one month, two month lead times
% essentially making the same prediction for every day in a month
CHad = NaN(length(lagRangeHad), 3*length(years_test));
for ii = 1:length(lagRangeHad)

	% hot days at specified lead time
	allMoYrLag = datenum(year(alldays), month(alldays) - lagRangeHad(ii) , 15);

	% index of hadsst data associated with that lead time
	loc = NaN(1,length(allMoYrLag));
	for ct = 1:length(allMoYrLag)
		loc(ct) = find(hadsst_comp.time == allMoYrLag(ct));
	end

	hadComp = nandetrend(nanmean(Zmat(:, loc), 2),'constant');
	% project onto the anomalies at the same lead time
	idxPredict = find(ismember(year(hadsst_comp2.time), years_test) & ...
			ismember(month(hadsst_comp2.time), [6:8]-lagRangeHad(ii)));
	sstAnom = nandetrend(Zmat2(:,idxPredict),'constant');
	numvals = sum(~isnan(sstAnom), 1);
	sstAnom(isnan(sstAnom)) = 0; hadComp(isnan(hadComp)) = 0;
	CHad(ii, :) = hadComp'*sstAnom./(numvals - 1);

end

CHad = bsxfun(@times, CHad, 1./std(CHad,[],2));

% transform CHad to have daily resolution for use in ROC scores
summerDaysPredict = bsxfun(@plus, datenum(years_test,1,1), summer_doy' - 1);
summerDaysPredict = summerDaysPredict(:);
summerMoPredict = unique([year(summerDaysPredict), month(summerDaysPredict)], 'rows');
summerMoPredict = [summerMoPredict(:,2) summerMoPredict(:,1)];
CHadDaily = NaN(length(lagRangeHad), length(summerDaysPredict));
for ct = 1:length(summerMoPredict) % match days to monthly predictions
	idxMo = find(ismember([month(summerDaysPredict), year(summerDaysPredict)], ...
		summerMoPredict(ct,:), 'rows'));
	CHadDaily(:, idxMo) = repmat(CHad(:,ct),[1 length(idxMo)]);
end

lagsToPlot = [2]; % plot 2 months
[roc_scoreHad,TPRHad, FPRHad, cutoffsHad] = pepROC(lagRangeHad, lagsToPlot, CHadDaily, 'HadPEP', ...
	summer_doy, alldaysHad, predictandName, ...
	GHCND_clustered_had, fig_folder);

%% And calculate null bounds for HadSST3 ROC scores
rng('default')
clear cache_input
summerTime = GHCND_clustered_had.time(ismember(GHCND_clustered_had.doy, summer_doy));
nsurr = 1e4;
rocScoreSurrHad = NaN(nsurr, length(lagRangeHad));
makeplots = 0;
cache_input.nsurr = nsurr;
cache_input.alldays = alldaysHad;
cache_input.lag_range = lagRangeHad;
cache_input.metric = 'alldaysHad';
cache_input.C = CHadDaily;
ROCsurrHash = DataHash(cache_input);
matname = [cacheDir '/ROCsurr' ROCsurrHash '.mat'];
if exist(matname,'file')
	load(matname)
else
	for ct = 1:nsurr
		disp([num2str(ct)])
		idxSurr = datasample(1:nyrsHad, nyrsHad);
		if strcmp(cache_input.metric,'alldaysHad')
			hotMat = ismember(GHCND_clustered_had.time, alldaysHad);
		elseif strcmp(cache_input.metric,'startDateHad')
			hotMat = ismember(GHCND_clustered_had.time, startDateHad);
		end
		hotMat = reshape(hotMat(ismember(GHCND_clustered_had.doy, summer_doy)),[60 nyrsHad]);
		hotMatSurr = hotMat(:, idxSurr);
		predictandSurr = summerTime(hotMatSurr(:));
		[rocScoreSurrHad(ct, :)] = ...
			calcROC(GHCND_clustered_had, CHadDaily, [], lagRangeHad, predictandSurr, ...
			summer_doy, makeplots, fig_folder, 'pacific', cutoffPerc);

	end
	save(matname,'rocScoreSurrHad');
end

% calculate p-value for ROCs
for ct = 1:length(lagRangeHad)
	pvalHadSST(ct) = length(find(rocScoreSurrHad(:,ct) > roc_scoreHad(ct)))/nsurr;
end
