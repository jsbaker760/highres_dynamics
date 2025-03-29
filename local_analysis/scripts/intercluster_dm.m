function ClusterDM=intercluster_dm(DM,clusters,C)
% DM: a distance matrix of all isolates
% clusters: the number of the cluster to which each isolate belongs 
% C: The clusters to include in the calculation. Elements not present in "clusters" are excluded

% number of clusters
N = numel(C);

% initialize
ClusterDM=zeros(N);

% remove clusters not being considered

% boolean
keep = ismember(clusters,C);

% resize
DM=DM(keep,keep);
clusters = clusters(keep);

% remove self-self comparisons in DM by turning them into NaN
DM(eye(size(DM))==1)=NaN;

%get inter-cluster distances

% arrays of logicals, true where isolate cluster is c in C
logical_clusters = arrayfun(@(c) clusters == c, C, 'UniformOutput', false);

% loop through each unique pair of lineages once and put values on both sides of matrix
for i = 1:N
    for j = i+1:N
        % Extract distances between clusters i and j
        distances = DM(logical_clusters{i}, logical_clusters{j});
        
        % Compute the median distance
        ClusterDM(i, j) = median(distances(:), 'omitnan');
        ClusterDM(j, i) = ClusterDM(i, j); % put values in both since symmetrical
    end
end


