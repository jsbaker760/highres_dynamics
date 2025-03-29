%% load data structures for clustering
% these are the isolate-wise alignment data for all samples aligned to
% Pacnes_C1 and SepepidermidisATCC12228
% Some variables were consolidated to save space, ex.
% coverage_forward_strand & coverage_forward_strand are summed to create coverage
Sepi_isolate_mutations = load('data/SepidermidisATCC12228_UnfilteredMutationVars.mat');
Cacnes_isolate_mutations = load('data/Pacnes_C1_UnfilteredMutationVars.mat');

% These booleans are true where an isolate is a C. acnes or S. epidermidis
% isolate, as determined by it's core genome content. 
load data/is_cacnes_isolate.mat
load data/is_sepi_isolate.mat

% These structures contain consolidated information for 
% all analyzed instances of clustering for both species and for
% both DBSCAN and ANCHOR
% they are generated with the script consolidate_cluster_data.m
Cacnes = load('CacnesClusterPlotting.mat');
Sepi = load('SepiClusterPlotting.mat');

% these are the distance matrices which were actually used for clustering
% and were generated using 
load data/distance_matrices.mat CacnesDM SepiDM

% These are the best instances of clustering for both species, 
% exactly as output by the algorithms
load data/clusters_raw.mat cacnes_clusters_raw sepi_clusters_raw

% Since cluster numbrs are arbirary, and clusters were numbered in a
% certain way previously, the sets of best clusters were re-numbered to
% match and prevent confusion
sepi_clusters_renumbered=load('sepi_clusters.mat','clusters');sepi_clusters_renumbered=sepi_clusters_renumbered.clusters;
cacnes_clusters_renumbered=load('cacnes_clusters.mat','clusters');cacnes_clusters_renumbered=cacnes_clusters_renumbered.clusters;

% *clusters_renumbered includes singletons and doubletons, but we decided
% to exclude these in practice.
cacnes_clusters = cacnes_clusters_renumbered;
sepi_clusters = sepi_clusters_renumbered;
cacnes_clusters(arrayfun(@(x) (x>0&(sum(cacnes_clusters==x)<=2)), cacnes_clusters))=0;
sepi_clusters(arrayfun(@(x) (x>0&(sum(sepi_clusters==x)<=2)), sepi_clusters))=0;

%% Calculate local clustering coefficients

cacnes_lcc_all=calculate_lcc(CacnesDM,10);
sepi_lcc_all=calculate_lcc(SepiDM,10);

% lineage-level data for within-cluster SNPs
load data/LineageData.mat

%% Histogram of pairwise distances
cacnes_samelin=(cacnes_clusters>0&cacnes_clusters'>0)&(tril(ones(size(CacnesDM)),-1)==1)&(cacnes_clusters==cacnes_clusters');
cacnes_difflin=(cacnes_clusters>0&cacnes_clusters'>0)&(tril(ones(size(CacnesDM)),-1)==1)&(cacnes_clusters~=cacnes_clusters');
sepi_samelin=(sepi_clusters>0&sepi_clusters'>0)&(tril(ones(size(SepiDM)),-1)==1)&(sepi_clusters==sepi_clusters');
sepi_diflin=(sepi_clusters>0&sepi_clusters'>0)&(tril(ones(size(SepiDM)),-1)==1)&(sepi_clusters~=sepi_clusters');

f=figure;
subplot(2,1,1)
bins=0:1000:(5*10^4);
bins_same =histcounts(CacnesDM(cacnes_samelin),bins);
bins_dif =histcounts(CacnesDM(cacnes_difflin),bins);
b = bar(bins(2:end),([bins_same;bins_dif])','stacked');
pbaspect([1 1 1])
xlim([0 5*10^4])
xlabel('pairwise SNP distances, C. acnes isolates')
ylabel('number of pairs')

subplot(2,1,2)
bins=0:1000:(10^5);
bins_same =histcounts(SepiDM(sepi_samelin),bins);
bins_dif =histcounts(SepiDM(sepi_diflin),bins);
b = bar(bins(2:end),([bins_same;bins_dif])','stacked');
pbaspect([1 1 1])
xlim([0 5*10^4])
xlabel('pairwise SNP distances, S. epidermidis isolates')
ylabel('number of pairs')

f.Renderer="painters";


%% LCC plots
% All isolates have lower connectivity at low E, but it's worse for bad
% isolates
% when you remove these bad isolates, it improives connectivity for
% clustered isolates
% however, they're still not perfect
f=figure;
subplot(2,1,1)
lcc = cacnes_lcc_all;
scatter(lcc.E_all_isolates(:),lcc.LCC_all_isolates(:),1,'filled','MarkerFaceColor',[0 0 0],'markerfacealpha',1)
xlabel('edge length (SNPs')
ylabel('Local clustering coefficient')
pbaspect([1 1 1])
subplot(2,1,2)
lcc = sepi_lcc_all;
scatter(lcc.E_all_isolates(:),lcc.LCC_all_isolates(:),1,'filled','MarkerFaceColor',[0 0 0],'markerfacealpha',1)
xlabel('edge length (SNPs')
ylabel('Local clustering coefficient')
pbaspect([1 1 1])

%% Reference genome positions in common amongst lineages
% The more distsntly related isolates are, the fewer genome positions they have in common. 
% this not only means that significantly less gene content can be compared amongst isolates, 
% it also means that there are non-linear differences in information
% density which distort distance metrics and make it difficult to cluster
% isolates.
f = figure;
subplot(1,2,1)
plot_refpositions_shared(Cacnes_isolate_mutations,cacnes_clusters,is_cacnes_isolate)
xlabel('number of C. acnes lineages being compared')
ylabel('variant positions comparable between lineages')
pbaspect([1 1 1])
xlim([0 100])
ylim([0 3*(10^5)])
subplot(1,2,2)
plot_refpositions_shared(Sepi_isolate_mutations,sepi_clusters,GoodStaphisolate'==1)
xlabel('number of S. epidermidis lineages being compared')
ylabel('variant positions comparable between lineages')
pbaspect([1 1 1])
xlim([0 100])
ylim([0 3*(10^5)])
ylim([10^5 3*10^5]) 
% methods: 100 subsets, each pair independently calculates positions not N
% between either


%% F5: pairwise genome distance vs core genome differences
% When comparing individual paira of lineage assemblies, more closely related lineages
% have more gene content in common. 
% the meaning is that there is less information to cluster when comparing
% distanly related isolates
ucacnes_clusters = unique(cacnes_clusters(cacnes_clusters>0));
usepi_clusters = unique(sepi_clusters(sepi_clusters>0));

CacnesClusterDM=intercluster_dm(CacnesDM,cacnes_clusters,ucacnes_clusters);
SepiClusterDM=intercluster_dm(SepiDM,sepi_clusters,usepi_clusters);

% lineage 105 excluded?
f=figure; hold on
subplot(2,1,1)
proportion_core_genes_shared(ucacnes_clusters,CacnesClusterDM,"cacnes");
xlabel('pairwise C. acnes core genome SNP distance');
ylabel('number of genes in common')
pbaspect([1 1 1])
subplot(2,1,2)
proportion_core_genes_shared(usepi_clusters,SepiClusterDM,"sepi");
xlabel('pairwise S. epidermidis core genome SNP distance');
ylabel('number of genes in common')
pbaspect([1 1 1])

%% Proportion of N SNPs (pairwise measurement) to N SNPs (core genes only)
% when simulating how many SNPs would be observed between two pairs of
% isolates alone in their shared gene content (versus the number of
% SNPs observed when comparing positions across isolates) there is a 
% multi-fold increase in actual-versus-apparent distances in closely
% related lineages, implicating that SNPs in non-core genes are 
% critical for differentiating lineages yet are missed when comparing all
% isolates.

f1 = figure;
snp_distance_distortion(Cacnes_isolate_mutations,cacnes_clusters,is_cacnes_isolate)
f2 = figure;
snp_distance_distortion(Sepi_isolate_mutations,sepi_clusters,is_sepi_isolate)


%% Histograms of allele frequencies 
% most of the samples which were not clustered
% were not obviously bad based on reference genome position
% so a filter would not be easy to apply
load('data/cdhit_filter.mat', 'CacnesCDHITfilter')
bool_all=CacnesCDHITfilter;
bool_kept=(is_cacnes_isolate>0);
f1=figure;
maf_example(Cacnes_isolate_mutations,cacnes_clusters_raw,bool_all,bool_kept)
ylim([0 .8])
load('data/cdhit_filter.mat', 'SepiCDHITfilter')
bool_all=SepiCDHITfilter;
bool_kept=bool_all;
f2=figure;
maf_example(Sepi_isolate_mutations,sepi_clusters_raw,bool_all,bool_kept)
ylim([0 .8])


%% Graph representation of excluding isolates (example from C. acnes)
% show that unclustered isolates have bad connectivity
% and their exclusion allows better clustering
show_graph_traversal(CacnesDM,cacnes_clusters_raw)


%% Get basic data for best cluster instances used
idx_CacnesANCHOR_best = Cacnes.n_clustered_new==max(Cacnes.n_clustered_new);
idx_CacnesDBSCAN_best = Cacnes.is_anchor==0&Cacnes.n_clustered_new==max(Cacnes.n_clustered_new(Cacnes.is_anchor==0));
idx_SepiANCHOR_best = Sepi.is_anchor==1&Sepi.n_clustered_new==max(Sepi.n_clustered_new(Sepi.is_anchor==1));
idx_SepiDBSCAN_best = Sepi.is_anchor==0&Sepi.n_clustered_new==max(Sepi.n_clustered_new(Sepi.is_anchor==0));

cacnes_match = sum(Cacnes.clusters_new_all==cacnes_clusters_raw,2)==numel(cacnes_clusters_raw);
sepi_match = sum(Sepi.clusters_new_all==sepi_clusters_raw,2)==numel(sepi_clusters_raw);

Cacnes.instance_best_algo = idx_CacnesDBSCAN_best|idx_CacnesANCHOR_best;
Sepi.instance_best_algo = idx_SepiANCHOR_best|idx_SepiDBSCAN_best;
Cacnes.instance_used = cacnes_match;
Sepi.instance_used = sepi_match;

% for S. epidermidis, there are two best solutions which were found 3x each
sepi_ANCHORbest = Sepi.clusters_new_all(idx_SepiANCHOR_best,:);
sepi_ANCHORbest_u= unique(sepi_ANCHORbest,'rows');
sepi_DBSCANbest = Sepi.clusters_new_all(idx_SepiDBSCAN_best,:);
% DBSCAN has 48 unique optimal solutions, first one was chosen
% why? not super clear. all have same W
sepi_DBSCANbest_u= unique(sepi_DBSCANbest,'rows');
[NSepi,WSepi,CSepi] = getNWC(SepiDM,sepi_DBSCANbest_u);
% cacnes anchor solution (24) was exactly the same across all
cacnes_ANCHORbest = Cacnes.clusters_new_all(idx_CacnesANCHOR_best,:);
cacnes_ANCHORbest_u = unique(cacnes_ANCHORbest,'rows');
% in DBSCAN, 81 exactly the same instances were found
cacnes_DBSCANbest = Cacnes.clusters_new_all(idx_CacnesDBSCAN_best,:);
cacnes_DBSCANbest_u = unique(cacnes_ANCHORbest,'rows');

%% margin vs N clustered, colored by cutoff coefficient percent
% show optimal instances of clustering
f1=plot_lcc_nclustered(Cacnes);
subplot(2,1,1)
ylim([1800 2250])
yticks([1800:50:2250])
subplot(2,1,2)
ylim([1800 2250])
yticks([1800:50:2250])
f2=plot_lcc_nclustered(Sepi);
subplot(2,1,1)
ylim([1700 2150])
yticks([1700:50:2150])
subplot(2,1,2)
ylim([1700 2150])
yticks([1700:50:2150])

%% Number of clustered isolates vs number of clusters
f1=nclusters_vs_nclustered(Cacnes);
subplot(2,1,1)
ylim([1800 2250])
xlim([90 140])
yticks([1800:50:2250])
subplot(2,1,2)
ylim([1800 2250])
xlim([90 140])
yticks([1800:50:2250])

f2=nclusters_vs_nclustered(Sepi);
subplot(2,1,1)
ylim([1700 2150])
xlim([100 200])
yticks([1700:50:2150])
subplot(2,1,2)
ylim([1700 2150])
xlim([100 200])
yticks([1700:50:2150])

%% Plot clustering instances by parameters
% and show that greedy addition improves things
f=plot_clusterparams(Cacnes)
subplot(2,2,1);
ylim([1800 2250])
subplot(2,2,2);
ylim([1800 2250])
subplot(2,2,3);
zlim([1800 2250])
subplot(2,2,4);
zlim([1800 2250])


f=plot_clusterparams(Sepi)
subplot(2,2,1);
ylim([1600 2150])
subplot(2,2,2);
ylim([1600 2150])
subplot(2,2,3);
zlim([1600 2150])
subplot(2,2,4);
zlim([1600 2150])

f1=plot_greedy_addition_results(Cacnes)
subplot(2,1,1);
ylim([1800 2250])
xlim([1800 2250])
subplot(2,1,2);
ylim([1800 2250])
xlim([1800 2250])
 
f2=plot_greedy_addition_results(Sepi)
subplot(2,1,1);
ylim([1600 2150])
xlim([1600 2150])
subplot(2,1,2);
ylim([1600 2150])
xlim([1600 2150])
