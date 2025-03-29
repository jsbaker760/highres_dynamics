%% This script consolidates clustering data and does the "greedy addition" step

%% initialize
% this script was called in an HPC environment for each species
% with different input variables
% but this logic statement has been added after-the-fact
% to make one wrapper for both species
for species = 1:2
    if species==1
        load data/SepidermidisATCC12228_clustering_instances.mat
        Species = SepidermidisATCC12228;
        name="sepi";
    elseif species==2
        load data/Pacnes_C1_clustering_instances.mat
        Species = Pacnes_C1;
        name="cacnes";
    end
    % reconstruct indices used in filtering
    Ucoeficients= .5:.01:1;
    Umargins = [1:1:100 200:100:1000];
    coef_combinations = [numel(Ucoeficients) numel(Umargins)];
    C = prod(coef_combinations);
    [a,b]=ind2sub(coef_combinations,1:C);
    co = Ucoeficients(a)';
    marid = (Umargins(b)');
    % where unique rows have multiple values, put them into cells
    margins_all = marid;
    coefficients_all = co;
    margins_urows=arrayfun(@(x) {margins_all(Species.CoefRows.ic==x)}, 1:max(Species.CoefRows.ic));
    coefficients_urows=arrayfun(@(x) {coefficients_all(Species.CoefRows.ic==x)}, 1:max(Species.CoefRows.ic));
    n_included = sum(Species.CoefRows.UrowsMeetsCutoff,2);

    % get clusters from before greedy addition
    ANCHOR_before = Species.ANCHOR.clusters_unique_consistent;
    DBSCAN_before = Species.DBSCAN.clusters_unique_consistent;

    % get number of included isolates from each filter
    DB_n_included=n_included(Species.DBSCAN_coefrowidx);
    ANCHOR_n_included=n_included(Species.ANCHOR_coefrowidx);

    % get margins for each row
    margins_db = margins_urows(Species.DBSCAN_coefrowidx);
    margins_anchor = margins_urows(Species.ANCHOR_coefrowidx);

    % coefs for each row
    coefs_db = coefficients_urows(Species.DBSCAN_coefrowidx);
    coefs_anchor = coefficients_urows(Species.ANCHOR_coefrowidx);

    % resize above for consistency
    DB_n_included=DB_n_included(Species.DBSCAN.consistency);
    ANCHOR_n_included=ANCHOR_n_included(Species.ANCHOR.consistency);
    margins_db=margins_db(Species.DBSCAN.consistency);
    coefs_db=coefs_db(Species.DBSCAN.consistency);
    margins_anchor=margins_anchor(Species.ANCHOR.consistency);
    coefs_anchor=coefs_anchor(Species.ANCHOR.consistency);

    % algorithm specific params
    % ANCHOR
    C1 = [Species.C1(Species.ANCHOR.consistency)];

    % DBSCAN
    epsilon = [Species.epsilon(Species.DBSCAN.consistency)];
    minpts = [Species.minpts(Species.DBSCAN.consistency)];

    % boolean for below, true where is anchor instance (rest are DBSCAN)
    is_anchor = [ones(size(ANCHOR_before,1),1) ; zeros(size(DBSCAN_before,1),1)];

    % put them all together so in same size as clusters_new
    clusters_before =       [ANCHOR_before ; DBSCAN_before];
    n_included_both =       [ANCHOR_n_included; DB_n_included];
    margins_filter_both =   [margins_anchor'; margins_db'];
    coefs_filter_both =     [coefs_anchor'; coefs_db'];
    n_not_excluded =        [ANCHOR_n_included ; DB_n_included];
    C1 =                    [C1 ;  nan(size(DBSCAN_before,1),1)];
    epsilon =               [nan(size(ANCHOR_before,1),1) ; epsilon];
    minpts =                [nan(size(ANCHOR_before,1),1) ; minpts];

    % Too slow to do clustering for all instances, so subsample them
    if name=="cacnes"
        bool_exclude_rows = (sum(clusters_before>0,2)<1800)|max(clusters_before,[],2)<96|max(clusters_before,[],2)>150;
    elseif name=="sepi"
        bool_exclude_rows = sum(clusters_before>0,2)<1600|max(clusters_before,[],2)<100|max(clusters_before,[],2)>200;
    end

    % downsample data
    clusters_before(bool_exclude_rows,:)=[];
    is_anchor(bool_exclude_rows)=[];
    C1(bool_exclude_rows)=[];
    epsilon(bool_exclude_rows)=[];
    minpts(bool_exclude_rows)=[];
    n_not_excluded(bool_exclude_rows)=[];
    coefs_filter_both(bool_exclude_rows)=[];
    margins_filter_both(bool_exclude_rows)=[];
    n_included_both(bool_exclude_rows)=[];

    % check for consistency
    if name=="cacnes"
        CacnesClustersNew = append_clusters(CacnesDM,clusters_both);
        clusters_new_all=CacnesClustersNew;
        is_consistent = arrayfun(@(x) check_consistency(CacnesDM,clusters_new_all(x,:)), 1:size(clusters_new_all,1));
    elseif name=="sepi"
        SepiClustersNew = append_clusters(SepiDM,clusters_before);
        clusters_new_all=SepiClustersNew;
        is_consistent = arrayfun(@(x) check_consistency(SepiDM,clusters_new_all(x,:)), 1:size(clusters_new_all,1));

    end

    % remove inconsistent data
    clusters_new_all = clusters_new_all(is_consistent,:);
    clusters_before = clusters_before(is_consistent,:);
    n_not_excluded=n_not_excluded(is_consistent);
    n_included_both=n_included_both(is_consistent);
    minpts=minpts(is_consistent);
    margins_filter_both=margins_filter_both(is_consistent);
    is_anchor=is_anchor(is_consistent);
    epsilon=epsilon(is_consistent);
    coefs_filter_both=coefs_filter_both(is_consistent);
    C1=C1(is_consistent);
    %
    n_clustered_old=sum(clusters_before>0,2);
    n_clustered_new=sum(clusters_new_all>0,2);
    %
    if name=="cacnes"
        save('CacnesClusterPlotting.mat','C1','coefs_filter_both','epsilon','is_anchor','margins_filter_both','minpts','n_clustered_new','n_clustered_old','n_excluded','clusters_new_all','clusters_before','-v7.3');
    elseif name=="sepi"
        save('SepiClusterPlotting.mat','C1','coefs_filter_both','epsilon','is_anchor','margins_filter_both','minpts','n_clustered_new','n_clustered_old','n_excluded','clusters_new_all','clusters_before','-v7.3');
    end
end
%% helper to actually runs addition step
function newclusters = append_clusters(DM,clusters_old)

[U,S] = size(clusters_old);
newclusters = zeros(U,S);
clusters_old=parallel.pool.Constant(clusters_old);
DM = parallel.pool.Constant(DM);

parfor u=1:U
    n= add_back(DM.Value,clusters_old.Value(u,:),75);
    newclusters(u,:)=n;
end
end

