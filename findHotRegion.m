function [T95prime, T95] = findHotRegion(GHCND_clustered, baseline_percentile, ...
	station_areas, summerDays)
  
%% Get spatial 95th percentile
for i = 1:size(GHCND_clustered.anomTotal, 2)
	% sort stations by temperature
	[T index] = sort(GHCND_clustered.anomTotal(:, i), 'descend');
	% put station areas in same order
	A = station_areas(index);
	loc = find(~isnan(T));
	% find when the cumulative sum of areas is 5% of the domain
	areafrac = cumsum(A(loc))/sum(A(loc));
	[~,top5] = min(abs(areafrac - (100-baseline_percentile)/100));
	T95(i) = T(loc(top5));
end
T95summerAvg = mean(T95(summerDays));
T95prime = T95 - T95summerAvg;