function z_comp = remove_climatology(z_comp, varname)
% Remove climatology via 3 harmonics

nlat = numel(z_comp.lat);
nlon = numel(z_comp.lon);

dtime = z_comp.time(2) - z_comp.time(1);
if dtime == 1
	doy = z_comp.time - datenum(year(z_comp.time), 1, 1) + 1;
else
	disp('Need to fix HadSST3 code -- do it now!')
	keyboard
	doy = month(z_comp.time); % take monthly averages
end
 
if strcmp(varname,'sst')
	% switch missing values to NaN
	z_comp.data(z_comp.data < -10) = NaN;
end

% find leap days
leapDays = find(doy == 366);
ntime = numel(z_comp.time);
unique_doy = unique(doy);
timeTemp = z_comp.time;
timeTemp(leapDays) = [];
dataTemp = z_comp.data;
dataTemp(:,:,leapDays) = [];

% z_dev is deviation from climatology using three sinusoids
tBasis = ([1:365] - 0.5)/365;
basis1 = exp(2*pi*i*tBasis);
basis2 = exp(4*pi*i*tBasis);
basis3 = exp(6*pi*i*tBasis);

% estimate annual cycle using complete years
idxComplete = find(month(timeTemp) == 1 & day(timeTemp) == 1, 1, 'first'):...
	find(month(timeTemp) == 12 & day(timeTemp) == 31, 1, 'last');
loc = ismember(1:length(timeTemp), idxComplete);
zComplete = dataTemp;
zComplete(:,:,~loc) = []; % faster than direct subsetting
zComplete = reshape(zComplete,[nlon*nlat 365 length(idxComplete)/365]);

zClim = bsxfun(@minus, zComplete, mean(zComplete, 2)); % remove mean of each year
zClim = nanmean(zClim, 3); % nlon*nlat x 365

eof1 = 2/length(tBasis)*sum(bsxfun(@times, zClim, basis1), 2);
eof2 = 2/length(tBasis)*sum(bsxfun(@times, zClim, basis2), 2);
eof3 = 2/length(tBasis)*sum(bsxfun(@times, zClim, basis3), 2);

rec = real(bsxfun(@times, conj(eof1), basis1)) + ...
	real(bsxfun(@times, conj(eof2), basis2)) + ...
	real(bsxfun(@times, conj(eof3), basis3));

% remove climatology
data = reshape(z_comp.data,[nlon*nlat ntime]);
z_dev = NaN(size(data));
for ct = 1:length(doy)
	if doy(ct) ~= 366
		z_dev(:,ct) = data(:,ct) - rec(:,doy(ct));
	else
		z_dev(:,ct) = data(:,ct) - rec(:,365);
	end
end

if strcmp(varname, 'sst')
    % remove great lakes - strong relationship to east coast land temperature
    [LAT LON] = meshgrid(z_comp.lat,z_comp.lon);
    great_lakes = LON >268 & LON< 284 & LAT > 40 & LAT < 49;
    z_dev(great_lakes(:), :) = NaN;
end

z_dev = reshape(z_dev, size(z_comp.data));
% add zero mean for domain of interest (?)
z_comp.anom = z_dev;