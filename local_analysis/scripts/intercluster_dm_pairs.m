function ClusterDM=intercluster_dm_pairs(DM,clusters,C)

%% Initialize
N = numel(C);
ClusterDM=zeros(N);

%% remove clusters not being considered

keep = ismember(clusters,C);
DM=DM(keep,keep);
clusters = clusters(keep);

%% for each pair of lineages, get median distance

for i=1:N
    for j=1:N
        if i~=j
            ClusterDM(i,j)=max(max([DM(clusters==C(i),clusters==C(j))]));
        end
    end
end
