%% These two functions are the same, but they run clustering for each species seperately
% to save time. They are both run (on a Slurm based HPC) with the command
% run.sh

cd scripts/clustering/

% the input is DMs.mat, in the same folder as these two scripts and run.sh
SepidermidisATCC12228_ParameterSweep.m
Pacnes_C1_ParameterSweep.m
parameter_sweep_best_dm.m
% The outputs are Pacnes_C1.mat and SepidermidisATCC12228.mat


%% This is the greedy addition step. 


% The inputs are Pacnes_C1.mat and SepidermidisATCC12228.mat
greedy_addition
% The outputs are CacnesClustersNew.mat and SepiClustersNew.mat


%% Gets the best clustering instance from each


[Ncacnes,Wcacnes,Ccacnes] = getNWC(CacnesDM,CacnesClustersNew);
[Nsepi,Wsepi,Csepi] = getNWC(CacnesDM,CacnesClustersNew);

%%


%PREFILTERING
Ucoeficients= .5:.01:1;
Umargins = [1:100 200 500 1000];
midx = 1:numel(Umargins);
M = numel(midx);
coef_combinations = [numel(Ucoeficients) numel(Umargins)];
C = prod(coef_combinations);
[a,b]=ind2sub(coef_combinations,1:C);
co = Ucoeficients(a)';
marid = (midx(b)');

%
load GoodDMs.mat

% find the prefilters used
coefs=zeros(C,size(CacnesDM,1));
for i = 1:size(coefs,1)
    c= ACcoef(CacnesDM,marid(i));
    coefs(i,:)=c>=co(i);
    i/size(coefs,1)
end
%
for i = 1:numel(best_rows_original)
    coefidx=Pacnes_C1.ANCHOR_coefrowidx(best_rows_original(i))
    cutoff_used = Pacnes_C1.CoefRows.UrowsMeetsCutoff(coefidx,:);
    same_idx=find(mean(coefs==cutoff_used,2)==1)
    co(same_idx)
    marid(same_idx)
    foo=1
end

%% Get numbers of clusters


% IMPORTANT: the only reason this was done is because some other projects used
% lineage numbers for other projects, so this re-numbers
% the new clusters to match the new ones (it doesn't change cluster membership, of course)
oldclusters = load('data/CacnesClustersCurrent.mat');
oldclusters=oldclusters.clusters;
newclusters=load('CacnesClustersBestNew.mat');
newclusters=newclusters.clusters;

[clusters] =  getNewClusterIDs(oldclusters,newclusters);
save('data/cacnes_clusters.mat','clusters')

% For s. epi there's a different number of samples, so you need to
% corrospond them 
oldclusters = load('data/SepiClustersCurrent.mat');
oldclusters=oldclusters.clusters;
newclusters=load('SepiClustersBestNew.mat');
newclusters=newclusters.clusters    ;
b=load('GoodSepiBool23172.mat');
b=b.bool;
resized = zeros(size(b));
resized(b)=oldclusters;
oldclusters=resized;
[clusters] =  getNewClusterIDs(oldclusters,newclusters);
save('data/sepi_clusters.mat','clusters')


%%


SepiLineageIndices = load('SepiClustersBestNew.mat');SepiLineageIndices=SepiLineageIndices.clusters;
CacnesLineageIndices = load('CacnesClustersBestNew.mat');CacnesLineageIndices=CacnesLineageIndices.clusters;
CacnesCDHITfilter = load('data/CacnesCDHITfilter.mat');CacnesCDHITfilter=CacnesCDHITfilter.CDHITfilter;
SepiCDHITfilter = load('data/SepiCDHITfilter.mat');SepiCDHITfilter=SepiCDHITfilter.CDHITfilter;
SamplesNamesAll = load ('data/SampleNamesAll.mat');SamplesNamesAll=SamplesNamesAll.SampleNames;
DMs = load('data/DMs.mat');DMs=DMs.DMs;
[CacnesLineageIndices, CacnesIsClustered]=resize_clusters(CacnesLineageIndices,CacnesCDHITfilter,DMs.Pacnes_C1.excludedSamplesVariable{2});
[SepiLineageIndices, SepiIsClustered]=resize_clusters(SepiLineageIndices,SepiCDHITfilter,DMs.SepidermidisATCC12228.excludedSamplesVariable{3});
cladestring=getcladestring(CacnesLineageIndices,CacnesIsClustered,SepiLineageIndices,SepiIsClustered);

%%

plot_clustering_qc.m