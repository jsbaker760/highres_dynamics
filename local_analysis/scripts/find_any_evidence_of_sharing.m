function [CacnesEither,SepiEither]=find_any_evidence_of_sharing(IT,M,cutoff_isolates,cutoff_assigned,cutoff_real)
%% Get abundances to use in comparison for both species
plothandles = struct;
MetagenomicsCacnesAbundance = M.CombinedCacnesLineages;
MetagenomicsSepiAbundance=M.CombinedSepiLineages;

%% get cluster numbers
MetagenomicsClustersCacnes=M.CacnesLineageNumbers;
MetagenomicsClustersSepi=M.SepiLineageNumbers;

%% comparing subject-time (samplings), so get this for each row in M
MetagenomicsSubjectTime=arrayfun(@(x) strjoin([M.SID(x) string(M.TP(x))],'') , 1:size(M,1))';

%% run function to get numbers
CacnesIsolateNumber=mg2isolates(MetagenomicsClustersCacnes,MetagenomicsSubjectTime,IT(IT.SpeciesString=="cacnes",:));
SepiIsolateNumber=mg2isolates(MetagenomicsClustersSepi,MetagenomicsSubjectTime,IT(IT.SpeciesString=="sepi",:));

%% Find isolates found in EITHER data source (union)

CacnesIsolateNumber(sum(CacnesIsolateNumber,2)<cutoff_isolates,:)=0;
SepiIsolateNumber(sum(SepiIsolateNumber,2)<cutoff_isolates,:)=0;

MetagenomicsCacnesAbundance(sum(MetagenomicsCacnesAbundance,2)<cutoff_assigned,:)=0;
MetagenomicsSepiAbundance(sum(MetagenomicsSepiAbundance,2)<cutoff_assigned,:)=0;

MetagenomicsCacnesAbundance(MetagenomicsCacnesAbundance<cutoff_real)=0;
MetagenomicsSepiAbundance(MetagenomicsSepiAbundance<cutoff_real)=0;

CacnesEither=CacnesIsolateNumber>0|MetagenomicsCacnesAbundance>0;
SepiEither=SepiIsolateNumber>0|MetagenomicsSepiAbundance>0;

end

%% Get Isolate numbers array corrosponding to metagenomics abundances

function IsolateNumber=mg2isolates(Clusters,SubjectTime,IT)

% size of array to be made (isolate output same size as one from MG)
[I,J]=size(Clusters);
N=I*J;

SubjectTime=repmat(SubjectTime,1,J);

IsolateNumber=arrayfun(@(x) sum(IT.SubjectTime==SubjectTime(x)&IT.ClusterString==Clusters(x)), 1:N);

IsolateNumber=reshape(IsolateNumber,I,J);

end
