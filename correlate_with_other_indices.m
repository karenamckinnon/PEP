function [rPEP, r_predict, names] = correlate_with_other_indices(GHCND_clustered, C, summer_doy, yrs, nyrs, ...
	startDate, lag_range, T95prime, makeplots)
	
 
%% Interannual skill
CAvg = reshape(mean(C(lag_range >= 15 & lag_range <= 45,:), 1),[length(summer_doy) nyrs]);
CYr = mean(CAvg, 1);
% compare with median T95
T95summer = median(reshape(T95prime(ismember(GHCND_clustered.doy, summer_doy)),[length(summer_doy) nyrs]), 1);


rPEP = xcPH(CYr, T95summer, 1);

%% Correlate with other indices
indices_folder = '/n/huybers_lab/common/data/indices';
names = {'AO_monthly','PDO_monthly','NAO_monthly','PNA_monthly',...
	'ENSO_MEI','ENSO_monthly'};
summer_yrs = unique(year(GHCND_clustered.time(ismember(GHCND_clustered.doy, summer_doy))));
nyrs = length(summer_yrs);

% try to predict with other indices
for counter = 1:length(names)
	load([indices_folder '/' names{counter} '.mat'])
	index = eval(names{counter});
	if min(size(index.data)) > 1, disp(['Multiple indices in one structure for ' names{counter} '']); end
	doy = index.time - datenum(1,1,year(index.time)) + 1;

	% check predictive skill for preceding months
	for mo = [1:12];
		yrsUse = summer_yrs;
		if ~strcmp(names{counter},'ENSO_monthly')
			ts = index.data(month(index.time) == mo, 1);
		else
			ts = index.data(month(index.time) == mo, 7);
		end
		ts = ts(ismember(year(index.time(month(index.time) == mo)), yrsUse));
		r_predict(counter, mo) = xcPH(ts, T95summer(ismember(yrsUse, year(index.time(month(index.time) == mo)))), 1);

	end

end

