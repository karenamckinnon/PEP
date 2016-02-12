function precip_comp_censored = process_precip(precip_comp_censored, precip_fname);

nlat = length(precip_comp_censored.lat);
nlon = length(precip_comp_censored.lon);
ntime = length(precip_comp_censored.time);
precip_comp_censored.data = double(precip_comp_censored.data);
precip_comp_censored.info = ncinfo(precip_fname);
fillval = precip_comp_censored.info.Variables(4).Attributes(4).Value;
precip_comp_censored.data(precip_comp_censored.data == fillval) = NaN;

Zmat = reshape(precip_comp_censored.data,[nlon*nlat ntime]);
 
% some gridpoints at the northern edge of the domain are missing additional points
% for some reason. also set to zero for now
conus = ~isnan(Zmat(:,1));
nanidx = isnan(Zmat(conus,:));
dummy = Zmat(conus,:);
dummy(nanidx) = 0;
Zmat(conus, :) = dummy;

precip_comp_censored.data = reshape(Zmat,[nlon nlat ntime]);