function GHCND = get_ghcnd(yr_begin,yr_end, varname, ghcnd_matpath, station_latlims, station_lonlims)
% get_ghcnd(yr_begin,yr_end)
% Read in GHCND data for desired years and locations
%  
% Input: yr_begin = require data beginning on 1/1/yr_begin
%        yr_end = require that data goes through 12/31/yr_end
%        varname = 'TMAX' or 'TMIN'
%        ghncd_matpath = where are the GHCND matfiles?
%
% Output: GHCND = a structure with temperature data, location, station id,
% and time

% Check to see if the data has already been obtained
rng('default');
cache_inputs.yr_begin = yr_begin;
cache_inputs.yr_end = yr_end;
cache_inputs.varname = varname;
cache_inputs.ghcnd_matpath = ghcnd_matpath;
cache_inputs.station_latlims = station_latlims;
cache_inputs.station_lonlims = station_lonlims;
hash = DataHash(cache_inputs);
if exist(['/n/huybers_lab/common/ghcnd/analysis/composites/cache/ghcnddata' hash '.mat'],'file') % load existing?
    disp('All relevant data previously loaded... loading')
    load(['/n/huybers_lab/common/ghcnd/analysis/composites/cache/ghcnddata' hash '.mat'])
else % Read in matfiles
	% Require that data span summers on both start and end years
	t = datenum(yr_begin, 6, 1):datenum(yr_end, 8, 31);
	ntime = numel(t);
 
	% load metadata
	load([ghcnd_matpath 'metadata.mat'])
	fn = char(['metadata.' varname '_nobs']);
	fn1 = char(['metadata.' varname '_firstdate']);
	fn2 = char(['metadata.' varname '_lastdate']);
	fnlat = char('metadata.lat');
	fnlon = char('metadata.lon');
	date1 = datevec(eval(fn1));
	date2 = datevec(eval(fn2));

	% read in CONUS only
	conus_index =  eval(fnlat) > min(station_latlims) & eval(fnlat) < max(station_latlims) & ...
		eval(fnlon) > min(station_lonlims) & eval(fnlon) < max(station_lonlims);

	files_with_data = find(eval(fn)>0 & date1(:,1)' <= yr_begin & date2(:,1)' >= yr_end & conus_index);


	nf = length(files_with_data);
	location = NaN(nf, 3);
	data = NaN(nf, ntime);

	for ct = 1:nf
		fn = [ghcnd_matpath (metadata.station_codes{files_with_data(ct)})];
		SR = load(fn);
		stationRecord = SR.stationRecord;
		disp(['Station: ' stationRecord.meta.id ', ' num2str(ct) ' out of ' num2str(nf) ''])
		location(ct,:) = [metadata.lon(files_with_data(ct)) metadata.lat(files_with_data(ct)) metadata.elev(files_with_data(ct))];
		id(ct,:) = metadata.station_codes{files_with_data(ct)};
		% Get temperature and date vector. flag = 1 if fits date criteria

		[T d flag] = get_data(t,varname,stationRecord);

		% check again to make sure times match
		if flag == 1 && sum(abs(d - t)) ~=0, flag = 0; end;
		if flag
			data(ct,:) = T/10; % GHCND is in 0.1C integer units
		end
	end

	GHCND.location = location;
	GHCND.id = id;
	GHCND.data = data;
	GHCND.time = t;

    timestamp = datestr(now);
    save(['/n/huybers_lab/common/ghcnd/analysis/composites/cache/ghcnddata' hash '.mat'],'GHCND', 'timestamp');
end

return