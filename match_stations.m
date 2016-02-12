function [GHCND_clustered_had, match_idx] = match_stations(GHCND_QC_had, GHCND_clustered);
% Pick out the stations from the longer time period that are in the previous cluster

match_idx = ismember(GHCND_QC_had.id, GHCND_clustered.id, 'rows');
GHCND_clustered_had = subset(GHCND_QC_had, match_idx); 