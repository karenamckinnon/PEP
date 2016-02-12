function z_comp = removeTrend(z_comp);

if length(size(z_comp.data))~= 3, error('Wrong dimensions'); end

dataMat = reshape(z_comp.data,[size(z_comp.data,1)*size(z_comp.data,2) size(z_comp.data,3)]);
% DETREND removes the trend from each column of the matrix.
z_comp.anom = reshape(detrend(dataMat')',size(z_comp.data));

return 