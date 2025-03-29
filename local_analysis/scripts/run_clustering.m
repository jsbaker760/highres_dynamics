%% Move to directory with clustering scripts
% Methods: Clustering isolates into lineages

%% Run clustering for both species across a range of parameters
% This code calculates all instances first and then sub-samples them to 
% self-consistent instances so that they could be compared in QC steps
% and is therefore very slow.

% Parameters for each clustering algorithm are hoard-coded

% code requires distance_matrices.mat (input distance matrix for both species)

% designed to run in an HPC environment

% outputs (saves) Pacnes_C1_clustering_instances.mat and SepidermidisATCC12228_clustering_instances.mat,
% which are all the clustering instances for both species.
parameter_sweep_best_dm

%% Find isolates which were missed during clustering
% this step does the "greedy addition" in the Methods

% The inputs are Pacnes_C1_clustering_instances.mat and SepidermidisATCC12228_clustering_instances.mat,
% from parameter_sweep_best_dm.m

% designed to run in an HPC environment

% Outputs (saves) CacnesClusterPlotting.mat and SepiClusterPlotting.mat
% sor visualizing the results and searching for best instances
consolidate_clustering_data

%% get the best instances for both

Cacnes = load('data/CacnesClusterPlotting.mat');
Sepi = load('data/SepiClusterPlotting.mat');

CacnesMaxClustered = find(Cacnes.n_clustered_new==max(Cacnes.n_clustered_new),1);
SepiBestClustered = find(Sepi.n_clustered_new==max(Sepi.n_clustered_new),1);

cacnes_best = unique(Cacnes.clusters_new_all(CacnesMaxClustered,:),'rows');
sepi_best = unique(Sepi.clusters_new_all(SepiBestClustered,:),'rows');
% save these as data/cacnes_best.mat and data/sepi_best.mat 

%% make plots about clustering properties 

clustering_QC

%% Optional step: Renumber clusters indices
% %
% % IMPORTANT: 
% %
% % the purpose of this optional step is to re-number lineages
% % in this case, the arbitary labels of the clusters were re-numbered for
% % convenience.
% % This script does not change membership of clusters, and will break if 
% % there is not a perfect match in cluster membership
% oldclusters = load('data/relabel_cacnes_clusters.mat');
% oldclusters=oldclusters.clusters;
% newclusters=load('data/cacnes_best_instances.mat');
% newclusters=newclusters.clusters;
% [clusters] =  getNewClusterIDs(oldclusters,newclusters);
% save('data/cacnes_clusters.mat','clusters')
% 
% % For s. epi there's a different number of samples, so you need to
% % corrospond them 
% oldclusters = load('data/relabel_sepi_clusters.mat');
% oldclusters=oldclusters.clusters;
% newclusters=load('data/sepi_best_instances.mat');
% newclusters=newclusters.clusters;
% [clusters] =  getNewClusterIDs(oldclusters,newclusters);
% save('data/sepi_clusters.mat','clusters')
% 
% 
% %%
% 
