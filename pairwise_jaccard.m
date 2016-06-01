function D = pairwise_jaccard(XI,XJ)
% Pairwise jaccard distance, ignoring NaN entries.
% For the contingency table:
%    XI =|  1  |  0  |
% XJ = 1 |  a  |  c  |
% XJ = 0 |  b  |  d  |
%
% The Jaccard similarity is s = a/(a + c + b),
% with s on [0,1].
% I.e., it does not reward null associations.
% The Jaccard distance is d = 1 - s.
%
% This function simply computes the Jaccard distance
% in a pairwise fashion, ignoring NaNs.
% Input:  XI = 1xN vector (row of X for pdist)
%         XJ = mxN matrix (other rows of X for pdist)
% Output:  D = mx1 vector of distances.

% Strict method, which will fail (correctly) if inputs are not
% binary/logical:
%D = NaN(size(XJ,1),1);
%for j = 1:size(XJ,1)
%   ind = and(~isnan(XI),~isnan(XJ(j,:)));
%   D(j) = 1-sum(and(XI(ind),XJ(j,ind)))/sum(or(XI(ind),XJ(j,ind)));
%end

% ~30% faster method without input checking (input _must_ be binary with
% optional NaNs):
R0 = repmat(XI,size(XJ,1),1);
A = nansum(XJ.*R0,2);
D = 1 - A./(A+nansum(XJ.*abs(R0-1),2)+nansum(abs(XJ-1).*R0,2));
D(isnan(D)) = 1;
return