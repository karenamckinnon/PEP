function out = subset(in, index)
% For a given structure, pull out all data associated with a given index

out = in;
try out.location = out.location(index, :);end
try out.id = out.id(index, :); end
try out.data = out.data(index, :); end
try out.clusterid = out.clusterid(index); end
try out.anom = out.anom(index, :); end
try out.anomAnn = out.anomAnn(index, :); end
try out.anomTotal = out.anomTotal(index, :); end

return 