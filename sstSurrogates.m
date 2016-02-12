function [sstCompSurrBnds] = sstSurrogates(nsurr, sst_comp, summer_doy, alldays, ...
	summerLength, nyrs, prctilesSig, lag_range);

rng('default');
clear cache_input
summerTime = sst_comp.time(ismember(sst_comp.doy, summer_doy));
% summerSST = sst_comp.anom(:,:,ismember(sst_comp.doy, summer_doy));
hotMat = reshape(ismember(summerTime, alldays), [summerLength nyrs]);
cache_input.nsurr = nsurr;
cache_input.hotDays = alldays;
cache_input.prctilesSig = prctilesSig;
cache_input.lag_range = lag_range;
surrHash = DataHash(cache_input);
savename = ['/n/huybers_lab/common/ghcnd/analysis/composites/cache/sstSurr_' surrHash '.mat'];
if exist(savename,'file')
	disp('Loading SST surrogate bands')
	load(savename)
else 
	disp('Need to calculate SST surrogate bands')
	sstCompSurrBnds = NaN(size(sst_comp.anom,1), size(sst_comp.anom,2), length(prctilesSig), length(lag_range));
	for jj = 1:length(lag_range)
		disp([num2str(jj)])
		sstCompSurr = NaN(size(sst_comp.anom,1), size(sst_comp.anom,2), nsurr);
		lag = lag_range(jj);
		for kk = 1:nsurr

			surrInd = datasample(1:nyrs,nyrs);
			matSurr = hotMat(:, surrInd);

			dateUse = summerTime(matSurr(:));
			idxUse = find(ismember(sst_comp.time, dateUse)) - lag;
			idxUse(idxUse < 1) = [];
			idxUse(idxUse > size(sst_comp.anom,3)) = [];

			sstCompSurr(:,:,kk) = mean(sst_comp.anom(:,:,idxUse), 3);

		end
		sstCompSurrBnds(:, :, :, jj) = prctile(sstCompSurr, prctilesSig, 3);
	end


	timestamp = datestr(now);
	save(savename,'-v7.3','timestamp', 'sstCompSurrBnds','cache_input')
end

