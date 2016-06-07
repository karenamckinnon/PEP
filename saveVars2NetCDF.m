function saveVars2NetCDF(longfilename, dimnames, dims, varnames, vars)
% varnames and vars should be structures of the same size

% Open the file
ncid = netcdf.create(longfilename, 'NC_WRITE');
ndims = size(dimnames, 2);
nvars = size(varnames, 2);

% define dimensions
dimids = [];
for dimct = 1:ndims
	dimid{dimct} = netcdf.defDim(ncid, dimnames{dimct}, length(dims{dimct}));
	dimids = [dimids dimid{dimct}];
	dimvarid{dimct} = netcdf.defVar(ncid, dimnames{dimct}, 'double', [dimid{dimct}]);
end

for varct = 1:nvars
	varid{varct} = netcdf.defVar(ncid, varnames{varct}, 'double', dimids);
end

% We are done defining the NetCdf
netcdf.endDef(ncid);

% put dimensions in
for dimct = 1:ndims
	netcdf.putVar(ncid, [dimvarid{dimct}], dims{dimct});
end

% put variables in
for varct = 1:nvars
	netcdf.putVar(ncid, varid{varct}, vars{varct});
end

% close the netcdf
netcdf.close(ncid)
