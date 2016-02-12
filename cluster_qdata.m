function [clusters,c,Z] = cluster_qdata(T95,maxclust,dist_method,linkage_method,plot_dendro,plot_silhouette,figpath)
% Produces hierarchical agglomerative clusters of input data.
% Inputs: 
% T95 -- data; one row per station (for clustering on stations), with 0/1 timeseries. NaNs are
%        permitted.
% maxclust -- maximum number of clusters to return.
% dist_method -- 'jaccard' or 'dice'.
% linkage_method -- 'average', 'centroid', 'complete', 'median', 'single',
%                   'ward', or 'weighted'. See 'doc linkage' for more info.
%                   'average' is standard.
% plot_dendro -- 1 to plot a dendrogram
% plot_silhouette -- 1 to plot silhouette for each # of clusters (increases
% runtime by a factor of two!)
%
% Outputs:
% clusters -- columns of cluster IDs, with the i'th column representing
%             labels of the partitioning by i clusters.
% c        -- cophenetic coefficient (a measure of how well-clustered the
%             data are. 1 == perfect, 0 == terrible.
%
% The default call would be:
% [clusters,c] = cluster_qdata(T95,12,'jaccard','average',0,0)
X = T95;

if strcmpi(dist_method,'jaccard')
    Y = pdist(X,@pairwise_jaccard);
elseif strcmpi(dist_method,'dice')
    Y = pdist(X,@pairwise_dice);
end
Z = linkage(Y,linkage_method);

if plot_dendro
    figure('visible','off')
    dendrogram(Z,60)
    eval(['print -depsc ' figpath '/cluster_dendrogram.eps'])
end

c = cophenet(Z,Y);
clusters = cluster(Z,'maxclust',1:maxclust);

if plot_silhouette
    figure('visible','off')
    if strcmpi(dist_method,'jaccard')
        silhouette(X,clusters(:,maxclust),@pairwise_jaccard);
    elseif strcmpi(dist_method,'dice')
        silhouette(X,clusters(:,maxclust),@pairwise_dice);
    end
    drawnow
    eval(['print -depsc ' figpath '/cluster_silhouette_' num2str(maxclust) '.eps'])
end


return