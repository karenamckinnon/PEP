% Returns composite fields given a list of dates %
function z_comp = fetch_composite(ncfile,times_to_fetch,fldname,varargin)
 
z = netcdfobj(ncfile);
% Custom time conversion since NCEP has non-spec calendars with inaccurate
% metadata... otherwise, would use cdfdate2num
time_nc = z.vars.time.value;
time_units = z.vars.time.atts.units.value;
if strcmpi(time_units,'hours since 1-1-1 00:00:0.0')
    time = time_nc/24 - 2 + datenum('1-1-1 00:00:0.0');
elseif strcmpi(time_units,'hours since 1800-1-1 00:00:0.0')
	time = time_nc/24  + datenum('1800-1-1 00:00:0.0');
elseif strcmpi(time_units,'day') % case for OAFlux data
	time = datenum(1985,1,1):datenum(2009,12,31);
else
    time = cdfdate2num(time_units,'standard',time_nc);
end

% for monthly data, but all in terms of 15th of month
if strcmp(ncfile,'/n/huybers_lab/common/data/HadSST3/HadSST.3.1.1.0.median.nc')
    time = datenum(year(time_nc),month(time_nc),15);
elseif strcmp(ncfile,'/n/huybers_lab/common/data/HadSST2/5deg/HadSST2_1850on.nc')
    time = datenum(year(time),month(time),1);
end

% Get indices to load
ind = ismember(time, times_to_fetch);

% Load the field at those times
if nargin == 3
	data = ncread(ncfile, fldname);
else
	start = varargin{1};
	count = varargin{2}
	data = ncread(ncfile, fldname, start, count);
end
if any(size(data)==1), data = squeeze(data);end

z_comp.data = data(:,:,ind);
z_comp.lat = z.vars.lat.value;
z_comp.lon = z.vars.lon.value;
if isempty(z.vars.lat.value), z_comp.lat = z.vars.latitude.value; end
if isempty(z.vars.lon.value), z_comp.lon = z.vars.longitude.value; end

z_comp.time = time(ind);
z_comp.time_nc = time_nc(ind);

return