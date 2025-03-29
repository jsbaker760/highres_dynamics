function [T_all] = fill_outgroups(T,distance_matrices,CacnesUnfilteredMutationsFile,SepiUnfilteredMutationsFile)
%% Get names, clusters, and species of isolates for all clustered isolates

% load table with links to .fastqs
T = readtable(T,Delimiter=',');
% load distance matrices
load(distance_matrices,'CacnesDM','SepiDM');

% in lineages with more than 2 isolates, use coverage is for finding "best" isolate for outgroup
CacnesCov =load(CacnesUnfilteredMutationsFile,'coverage_unfiltered');
CacnesCov=mean(CacnesCov.coverage_unfiltered,1);

% get samplenames for both
SepiCov =load(SepiUnfilteredMutationsFile,'coverage_unfiltered','SampleNames_unfiltered');
AllSampleNames = string(SepiCov.SampleNames_unfiltered);
SepiCov=mean(SepiCov.coverage_unfiltered,1);

%% boolean for which rows are sepi and cacnes

is_cacnes = T.species=="cacnes";
is_sepi = T.species=="sepi";

%% resize

CacnesNames = T.SampleName(is_cacnes);
SepiNames = T.SampleName(is_sepi);

CacnesClusters = T.ClusterNumber(is_cacnes);
SepiClusters   = T.ClusterNumber(is_sepi);

%% pull out coverage

CacnesCov=arrayfun(@(x) CacnesCov(x==AllSampleNames), CacnesNames);
SepiCov=arrayfun(@(x) SepiCov(x==AllSampleNames), SepiNames);

%% get distance matrices to find outgroups
% note that since we removed singletons, we have to resize the distance
% matrix
idxCacnes = arrayfun(@(x) find(x==CacnesSampleNamesDM), CacnesNames);
idxSepi = arrayfun(@(x) find(x==SepiSampleNamesDM), SepiNames);

CacnesDM = CacnesDM(idxCacnes,idxCacnes); %#ok<*NODEF>
SepiDM = SepiDM(idxSepi,idxSepi);

%% find outgroup samples
outgroupclades_cacnes = get_outgroups(CacnesClusters,CacnesDM);
outgroupclades_sepi = get_outgroups(SepiClusters,SepiDM);

%% Get highest coverage representative
% first get all unique clade numbers
UCacnesClusters = unique(CacnesClusters);
uSepiClusters = unique(SepiClusters);

cluster_representatives_cacnes = CacnesNames(arrayfun(@(x) find(x==CacnesClusters&CacnesCov==(max(CacnesCov(x==CacnesClusters)))), UCacnesClusters));
cluster_representatives_sepi = SepiNames(arrayfun(@(x) find(x==SepiClusters&SepiCov==(max(SepiCov(x==SepiClusters)))), uSepiClusters));


%% add to table
TableCacnes = table;
idx=1;
for u = 1:numel(UCacnesClusters)
    ClusterNumber = UCacnesClusters(u);
    Members = CacnesNames(CacnesClusters==ClusterNumber);
    OutgroupIsolates = cluster_representatives_cacnes(outgroupclades_cacnes{u});
    for m = 1:numel(Members)
        Row = string(T.SampleName)==Members(m);
        TableCacnes.Path(idx)=T.Path(Row);
        TableCacnes.SampleName(idx)=T.SampleName(Row);
        TableCacnes.ReferenceGenome(idx)=strjoin(["cacnes" "clade" string(ClusterNumber)],'_');
        TableCacnes.ProviderNames(idx)=T.ProviderNames(Row);
        TableCacnes.is_outgroup(idx)=0;
        idx=idx+1;

    end
    for o = 1:numel(OutgroupIsolates)
        Row = string(T.SampleName)==OutgroupIsolates(o);
        TableCacnes.Path(idx)=T.Path(Row);
        TableCacnes.SampleName(idx)=T.SampleName(Row);
        TableCacnes.ReferenceGenome(idx)=strjoin(["cacnes" "clade" string(ClusterNumber)],'_');
        TableCacnes.ProviderNames(idx)=T.ProviderNames(Row);
        TableCacnes.is_outgroup(idx)=1;
        idx=idx+1;
    end
end

TableSepi = table;
idx=1;
for u = 1:numel(uSepiClusters)
    ClusterNumber = uSepiClusters(u);
    Members = SepiNames(SepiClusters==ClusterNumber);
    OutgroupIsolates = cluster_representatives_sepi(outgroupclades_sepi{u});
    for m = 1:numel(Members)
        Row = string(T.SampleName)==Members(m);
        TableSepi.Path(idx)=T.Path(Row);
        TableSepi.SampleName(idx)=T.SampleName(Row);
        TableSepi.ReferenceGenome(idx)=strjoin(["sepi" "clade" string(ClusterNumber)],'_');
        TableSepi.ProviderNames(idx)=T.ProviderNames(Row);
        TableSepi.is_outgroup(idx)=0;
        idx=idx+1;
    end
    for o = 1:numel(OutgroupIsolates)
        Row = string(T.SampleName)==OutgroupIsolates(o);
        TableSepi.Path(idx)=T.Path(Row);
        TableSepi.SampleName(idx)=T.SampleName(Row);
        TableSepi.ReferenceGenome(idx)=strjoin(["sepi" "clade" string(ClusterNumber)],'_');
        TableSepi.ProviderNames(idx)=T.ProviderNames(Row);
        TableSepi.is_outgroup(idx)=1;
        idx=idx+1;
    end
end
% add both to table
T_all = [TableCacnes ; TableSepi];
end

%% get outgroups for each cluster
function outgroups = get_outgroups(clusters,DM)

InterClusterDM = get_mean_intercluster_distances(clusters,DM);

% find 5 best
InterClusterDM(eye(size(InterClusterDM,1))==1)=NaN;
[~,idx] = sort(InterClusterDM,2,'ascend','MissingPlacement','last');
outgroups = arrayfun(@(x) {idx(x,1:5)}, 1:size(idx,1));

end

%
function InterClusterDM = get_mean_intercluster_distances(clusters,DM)

U = unique(clusters);

C = numel(U);

InterClusterDM = zeros(C);
for i = 1:(C)
    for j = 1:(C)
        InterClusterDM(i,j)=mean(DM(clusters==U(i),clusters==U(j)),'all');
    end
end

InterClusterDM(eye(size(InterClusterDM))==1)=0;
end