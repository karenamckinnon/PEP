function [clusters,cluster_counts] = relabel_clusters(clusters)
% Takes as input a single vector describing the cluster partitioning.
% Outputs a re-labeled clustering in which the clusters are labeled in
% descending order of number of members. Also returns a lookup table for
% the membership counts, original labels, and resulting labels.
 
% Sort clusters by size to use only the large ones
cluster_counts = NaN(max(clusters),2);
cluster_counts(:,1) = 1:max(clusters);
for j = 1:max(clusters)
   cluster_counts(j,2) = sum(clusters==j); 
end
cluster_counts = flipud(sortrows(cluster_counts,2));
% Relabel in descending order
cluster_counts = [cluster_counts (1:size(cluster_counts,1))'];
clusters_orig_labels = clusters;
for j = 1:size(cluster_counts,1)
   clusters(clusters_orig_labels==cluster_counts(j,1)) = cluster_counts(j,3); 
end
cluster_counts = cluster_counts(find(cluster_counts(:,2)>0),:);% keep only nonzero clusters (some initial labelings may screw this up otherwise)

return