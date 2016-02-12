% Download CPC data from online, and process into one netcdf file
% Note that, before 2006, the data is available in netcdf
% After 2006, each day must be downloaded as a binary file
 
years_to_download = 1979:2015;
loc_netcdf = 'ftp://ftp.cdc.noaa.gov/Datasets/cpc_us_precip';
dir_for_tmp = '/n/regal/huybers_lab/mckinnon/cpc-precip';

% years with netcdf
yrs = years_to_download(years_to_download <= 2006);

for ii = 1:length(yrs)
    fname = ['precip.V1.0.' num2str(yrs(ii)) '.nc'];
    if ~exist([dir_for_tmp '/' fname]),
    	system(['wget --directory-prefix=' dir_for_tmp '/ ' loc_netcdf '/' fname '']); % get files if not already there
    end
end

% get schema from last netcdf file
S = ncinfo([dir_for_tmp '/' fname]);

% lat/lon info from http://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/DOCU/PRCP_CU_GAUGE_V1.0CONUS_0.25deg.README
x = 230.125:0.25:(230.125+299*0.25);
y = 20.125:0.25:(20.125+119*0.25);
[X Y] = meshgrid(x,y);
X = X';
Y = Y';
npts = length(X(:));

% get other years
yrs2 = years_to_download(~ismember(years_to_download, yrs));
for ii = 1:length(yrs2)
    loc_bin = ['http://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/RT/' num2str(yrs2(ii)) ''];
    % get all files in directory
    dir2save = [dir_for_tmp '/' num2str(yrs2(ii)) ''];
    if ~exist(dir2save),
		mkdir(dir2save)
		system(['wget --directory-prefix=' dir2save '/ -r -l1 -nd --no-parent -A.RT.gz ' loc_bin '/']);
		system(['wget --directory-prefix=' dir2save '/ -r -l1 -nd --no-parent -A.RT ' loc_bin '/']);

		% go through each file and get data
		cd(dir2save)
		fnames = dir('PRCP_CU_GAUGE_V1.0CONUS_0.25deg.lnx.*.RT*');
		nfiles = length(fnames);
		% initialize
		precip = NaN(length(x), length(y), nfiles);
		time = NaN(nfiles,1);
		for jj = 1:nfiles
			fn = fnames(jj).name;
			if strcmp(fn(end-1:end), 'gz')
				system(['gunzip ' fn ''])
				fn2 = fn(1:end-3);
				% delete(fn) % automatically deletes .gz version
				fn = fn2;
			end
			time(jj) = datenum(str2num(fn(37:40)), str2num(fn(41:42)), str2num(fn(43:44)));
			f = fopen(fn);

			data = fread(f,[npts 2],'real*4','l'); % 4 byte little endian encoding
			vals = reshape(data(:,1),size(X));
			nstations = reshape(data(:,2),size(X));

			precip(:,:,jj) = vals/10; % original units are in 0.1 mm, want in mm

		end

		newnetcdfname = S.Filename;
		newnetcdfname(end-6:end-3) = num2str(yrs2(ii)); % replace year
		if exist(newnetcdfname), delete(newnetcdfname); end
		ncwriteschema(newnetcdfname,S)

		% switch time to same units as in netcdf file
		time_units = 'hours since 1800-01-01 00:00:0.0';
		time = (time - datenum(1800,1,1))*24; % hours since 1/1/1800

		% match missing value
		current_missing_val = min(precip(:));
		new_missing_val = S.Variables(4).Attributes(4).Value;
		precip(precip == current_missing_val) = new_missing_val;

		% write new data to file
		ncwrite(newnetcdfname,'time',time)
		ncwrite(newnetcdfname,'precip',precip)
    end
end

cd(dir_for_tmp)
% combine all netcdf files
save_dir = '/n/huybers_lab/common/data/cpc-precip';
if ~exist(save_dir), mkdir(save_dir); end
system(['ncrcat precip.V1.0.*.nc ' save_dir '/precip.V1.0.' num2str(min(years_to_download)) '.' num2str(max(years_to_download)) '.nc'])







