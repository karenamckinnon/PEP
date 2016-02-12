function GHCND_QC = remove_climo_ghcnd(GHCND_QC, years_to_average, cacheDir)
% years_to_average: 0 to remove linear trend, val > 0 to remove running mean

rng('default');
cache_inputs.id = GHCND_QC.id;
cache_inputs.years_to_average = years_to_average;
cache_inputs.time = GHCND_QC.time;
hash = DataHash(cache_inputs);

if exist([cacheDir '/anomdata' hash '.mat'],'file') % load existing?
    disp('Anomalies previously calculated')
    load([cacheDir '/anomdata' hash '.mat'])
    GHCND_QC.anomAnn = anomAnn;
	GHCND_QC.anomTotal = anomTotal;
else
	%% Remove climatology at each station using first three harmonics
	tBasis = ([1:365]-0.5)/365;
	basis1 = exp(2*pi*i*tBasis);
	basis2 = exp(4*pi*i*tBasis);
	basis3 = exp(6*pi*i*tBasis); 

	nstations = size(GHCND_QC.data, 1);
	ntime = length(GHCND_QC.time);

	leapDays = find(GHCND_QC.doy == 366);
	unique_doy = unique(GHCND_QC.doy);
	timeTemp = GHCND_QC.time;
	timeTemp(leapDays) = [];
	dataTemp = GHCND_QC.data;
	dataTemp(:,leapDays) = [];

	% get annual cycle using complete years
	idxComplete = find(month(timeTemp) == 1 & day(timeTemp) == 1, 1, 'first'):...
		find(month(timeTemp) == 12 & day(timeTemp) == 31, 1, 'last');
	nyrsUse = floor(length(idxComplete)/365);
	tsMat = reshape(dataTemp(:,idxComplete),[nstations 365 nyrsUse]);

	% remove annual mean
	tsClim = bsxfun(@minus, tsMat, nanmean(tsMat, 2));
	tsClim = nanmean(tsClim,3);
	if any(isnan(tsClim(:))), disp('NaNs in climatology'); keyboard; end

	eof1 = 2/length(tBasis)*sum(bsxfun(@times,tsClim,basis1),2);
	eof2 = 2/length(tBasis)*sum(bsxfun(@times,tsClim,basis2),2);
	eof3 = 2/length(tBasis)*sum(bsxfun(@times,tsClim,basis3),2);

	rec = real(bsxfun(@times, conj(eof1), basis1)) + ...
		real(bsxfun(@times, conj(eof2), basis2)) + ...
		real(bsxfun(@times, conj(eof3), basis3));

	% record max temperature day
	[~,maxTempDoy] = max(rec,[],2);

	recMat = NaN(size(GHCND_QC.data));
	for dayct = 1:(length(GHCND_QC.doy))
		if GHCND_QC.doy(dayct) ~= 366
			recMat(:, dayct) = rec(:, GHCND_QC.doy(dayct));
		else
			recMat(:, dayct) = rec(:, 365);
		end
	end

	anomAnn = GHCND_QC.data - recMat;

	if years_to_average == 0 % remove linear trend
		anomTotal = nandetrend(anomAnn')';
	else % remove annual running mean
		annmean = nanmoving_average(anomAnn, floor(365*years_to_average/2), 2);
		anomTotal = anomAnn - annmean;
	end

	timestamp = datestr(now);
	save([cacheDir '/anomdata' hash '.mat'],...
		'timestamp','anomAnn','anomTotal')
	GHCND_QC.anomAnn = anomAnn;
	GHCND_QC.anomTotal = anomTotal;
end
