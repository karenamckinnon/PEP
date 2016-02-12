function [startDate, endDate, totalTemp, alldays] = ...
	findHotDays(GHCND_clustered, cutofflength, T95prime, summerDays, ...
	yrs, summer_doy);
  
hot = T95prime > nanstd(T95prime(summerDays));

% first, get all dates
startDate = [];
endDate = [];
totalTemp = [];

nyrs = length(yrs);

for ct = 1:nyrs
	idxUse = ismember(year(GHCND_clustered.time), yrs(ct));
	timeUse = GHCND_clustered.time(idxUse);
	doyUse = GHCND_clustered.doy(idxUse);
	ts = hot(:,idxUse);
	tsSummer = ts(ismember(doyUse, summer_doy));
	timeSummer = timeUse(ismember(doyUse, summer_doy));
	dateStart = [];
	dateEnd = [];
	if tsSummer(1) == 1 && tsSummer(2) == 1,
		dateStart = timeSummer(1);
		indEnd = find(tsSummer == 0, 1, 'first');
		dateEnd = timeSummer(indEnd - 1);
	end
	for ii = 2:(length(tsSummer)-1)
		if tsSummer(ii-1) == 0 && tsSummer(ii) == 1 && tsSummer(ii+1) == 1
			dateStart = [dateStart;timeSummer(ii)];
			indEnd = find(tsSummer((ii+1):end) == 0, 1, 'first');
			if ~isempty(indEnd)
				dateEnd = [dateEnd; timeSummer(ii + indEnd - 1)];
			else
				dateEnd = [dateEnd; timeSummer(end)];
			end
		end
	end

	dtimeS = [NaN;dateStart(2:end) - dateEnd(1:end-1)];
	dtimeE = [dateStart(2:end) - dateEnd(1:end-1); NaN];

	% append events together
	dateStart(dtimeS < cutofflength) = [];
	dateEnd(dtimeE < cutofflength) = [];

	% Now, count up heat over time period
	totalTemp = NaN(length(dateStart),1);
	for jj = 1:length(dateStart)
		hotPeriod = ismember(timeSummer, dateStart(jj):dateEnd(jj));
		tempVals = T95prime(ismember(GHCND_clustered.time, timeSummer(hotPeriod)));
		tempVals(tempVals < 0) = 0;
		totalTemp(jj) = sum(tempVals);
	end

	startDate = [startDate;dateStart];
	endDate = [endDate;dateEnd];
	totalTemp = [totalTemp;totalTemp];
end

alldays = [];
for ct = 1:length(startDate)
	alldays = [alldays startDate(ct):endDate(ct)];
end
