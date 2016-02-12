function z_comp = removeClimoHadSST(z_comp, varname);

nlat = numel(z_comp.lat);
nlon = numel(z_comp.lon);
ntime = numel(z_comp.time);
dtime = z_comp.time(2) - z_comp.time(1);

% z_dev is deviation from climatology using three sinusoids
tBasis = ([1:12] - 0.5)/12;
basis1 = exp(2*pi*i*tBasis);
basis2 = exp(4*pi*i*tBasis);
basis3 = exp(6*pi*i*tBasis);
 
idxComplete = find(month(z_comp.time) == 1, 1, 'first'):...
	find(month(z_comp.time) == 12, 1, 'last');
loc = ismember(1:length(z_comp.time), idxComplete);
zComplete = z_comp.data;
zComplete(:,:,~loc) = []; % faster than direct subsetting
zComplete = reshape(zComplete,[nlon*nlat 12 length(idxComplete)/12]);

zClim = bsxfun(@minus, zComplete, nanmean(zComplete, 2)); % remove mean of each year
zClim = nanmean(zClim, 3); % nlon*nlat x 12

eof1 = 2/length(tBasis)*nansum(bsxfun(@times, zClim, basis1), 2);
eof2 = 2/length(tBasis)*nansum(bsxfun(@times, zClim, basis2), 2);
eof3 = 2/length(tBasis)*nansum(bsxfun(@times, zClim, basis3), 2);

rec = real(bsxfun(@times, conj(eof1), basis1)) + ...
	real(bsxfun(@times, conj(eof2), basis2)) + ...
	real(bsxfun(@times, conj(eof3), basis3));

% remove climatology
data = reshape(z_comp.data,[nlon*nlat ntime]);
mo = month(z_comp.time);
z_dev = NaN(size(data));
for ct = 1:length(mo)
	z_dev(:,ct) = data(:,ct) - rec(:,mo(ct));
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