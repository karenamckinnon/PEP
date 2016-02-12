function GHCND_QC = do_qc(GHCND)
% Remove stations that do not pass QC checks
% Check: at least 80% of summers must have at least 80% of their data
%    
% Input: GHCND = a structure with temperature data, location, station id,
% and time 
% Output: GHCND_QC of same format as GHCND, but with NaNs for stations that
% are not included

summer = find(month(GHCND.time) >= 6 & month(GHCND.time) <= 8);
summer_temperature = GHCND.data(:, summer);

nyrs = length(unique(year(GHCND.time)));
nsummerdays = length(summer)/nyrs;
nstations = size(summer_temperature, 1);

temp_mat = reshape(summer_temperature, [nstations nsummerdays nyrs]);
% fraction missing per summer (station x year)
frac_missing = squeeze(sum(isnan(temp_mat), 2))/nsummerdays;
% 80th percentile across years of the fraction of missing days
frac_missing_cutoff = prctile(frac_missing, 80, 2);
% only include stations where the 80th percentile is less than 0.2
% i.e. at least 80% of years must have 80% of data
idx_include = find(frac_missing_cutoff < 0.2);
GHCND_QC.data = GHCND.data(idx_include, :);
GHCND_QC.location = GHCND.location(idx_include, :);
GHCND_QC.id = GHCND.id(idx_include,:);
GHCND_QC.time = GHCND.time;

return
