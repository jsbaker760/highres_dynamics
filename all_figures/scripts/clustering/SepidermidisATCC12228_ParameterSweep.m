%% PARALLEL POOL
distcomp.feature( 'LocalUseMpiexec', false );
restoredefaultpath
parpool(12,'IdleTimeout',Inf);
dbstop if error

%%
load DMs.mat
DMs = DMs.SepidermidisATCC12228.DMsVariable;
D = numel(DMs);
%% PARAMETERS
%PREFILTERING
Ucoeficients= .5:.01:1;
Umargins = [1:100 200 500 1000];
% ANCHOR
Uneighborhoodr=[1:100];
% DBSCAN
uminpts=1:3;
uepsilon = [1:100];
minnum=1;
%% LOAD DATA
SepidermidisATCC12228=struct;
%%
midx = 1:numel(Umargins);
M = numel(midx);
coef_combinations = [numel(Ucoeficients) numel(Umargins)];
C = prod(coef_combinations);
[a,b]=ind2sub(coef_combinations,1:C);
co = Ucoeficients(a)';
marid = midx(b)';
%%
SepidermidisATCC12228.clustcoef = co;
SepidermidisATCC12228.cutoffcoef = marid;
%%
iterations = [numel(uepsilon) numel(uminpts)];
Udb = prod(iterations);
[e,m]=ind2sub(iterations,1:Udb);
epsilon =uepsilon(e);
minpts =uminpts(m);
%%
for d=1:D
    fprintf(['Starting sweep with DM' char(string(d)) ' of ' char(string(D)) '....\n'])
    %%
    distance_matrix = DMs{d};
    S = (size(distance_matrix,1));
    %% PREFILTERS
    % choose parameters
    coefs_all = zeros(M,S);
    dm=parallel.pool.Constant(distance_matrix);
    um = parallel.pool.Constant(Umargins);
    %% measuring coefs
    fprintf('measuring clustering coefficients....\n')
    parfor acc = 1:M
        coefs_all(acc,:) = ACcoef(dm.Value,um.Value(acc));
    end
    %%
    fprintf('choosing seed....\n')
    coef_margin_all = coefs_all(marid,:);
    meets_cutoff_all = coef_margin_all>= repmat(co,1,S);
    meets_cutoff_all(sum(meets_cutoff_all,2)<minnum,:)=[];
    seed_candidates = mean(coefs_all==1,1)==1;
    if sum(seed_candidates)>1
        seed_candidates=find(seed_candidates);
        seed_dm = sort(distance_matrix(seed_candidates,seed_candidates),2,'ascend');
        seed_candidates=seed_candidates(seed_dm(:,2)==max(seed_dm(:,2)));
        seed = false(S,1);seed(seed_candidates(1))=true;
    else
        seed = false(S,1);seed(seed_candidates)=true;
    end
    %%
    fprintf('getting unique rows....\n')
    [UrowsMeetsCutoff,ia,ic]=unique(meets_cutoff_all,'rows');
    CoefRows.UrowsMeetsCutoff=UrowsMeetsCutoff;
    CoefRows.ic=ic;
    CoefRows.ia=ia;
    CoefRows.seed=seed;
    fprintf('done with get_coef_rows....\n')
    %% RUNS ANCHOR
    fprintf('running Cacnes anchor....\n')
    %
    R = size(CoefRows.UrowsMeetsCutoff,1);
    Uanchor=numel(Uneighborhoodr);
    %
    [ANCHOR_clusters_all] = deal(zeros(R*Uanchor,S,'int16'));
    [ANCHORidx, C1_all]= deal(zeros(R*Uanchor,1));
    fprintf('starting ANCHOR iterations....\n')
    for r = 1:R
        idx2 = CoefRows.UrowsMeetsCutoff(r,:)==1;
        dm=parallel.pool.Constant(distance_matrix(idx2,idx2));
        seed=parallel.pool.Constant(CoefRows.seed(idx2));
        idx1 = (((r-1)*Uanchor)+1):(r*Uanchor);
        ANCHORidx(idx1)=r;
        C1_all(idx1)=Uneighborhoodr;
        I=numel(idx1);
        [rclusters_all]=deal(zeros(numel(idx1),sum(idx2),'int16'));
        Uneighborhood=parallel.pool.Constant(Uneighborhoodr);
        parfor i=1:I
            C1=Uneighborhood.Value(i);
            [anchors,clusters]=(AnchorClust(dm.Value,seed.Value,C1));
            rclusters_all(i,:)=clusters;
        end
        ANCHOR_clusters_all(idx1,idx2)=rclusters_all;
        fprintf([char(string((r/R)*100)) ' percent complete....\n'])
    end
    fprintf('finished ANCHOR iterations....\n')
    %% RUNS DBSCAN
    [DBSCAN_clusters_all]= deal(zeros(R*Udb,S,'int16'));
    [DBSCANidx, epsilon_all, minpts_all]=deal(zeros(R*Udb,1));
    fprintf('starting DBSCAN iterations....\n')
    for r = 1:R
        idx2 = CoefRows.UrowsMeetsCutoff(r,:)==1;
        dm=parallel.pool.Constant(distance_matrix(idx2,idx2));
        idx1 = (((r-1)*Udb)+1):(r*Udb);
        DBSCANidx(idx1)=r;
        epsilon_all(idx1)=epsilon;
        minpts_all(idx1)=minpts;
        I=numel(idx1);
        [rclusters_all]=deal(zeros(numel(idx1),sum(idx2),'int16'));
        Sr = parallel.pool.Constant(sum(idx2));
        parepsilon = parallel.pool.Constant(epsilon);
        parminpts = parallel.pool.Constant(minpts);
        parfor i = 1:I
            clusters = DBSCAN(dm.Value,parepsilon.Value(i),parminpts.Value(i));
            rclusters_all(i,:)=clusters;
        end
        DBSCAN_clusters_all(idx1,idx2)=rclusters_all;
        fprintf([char(string((r/R)*100)) ' percent complete....\n'])
    end
    fprintf('finishing DBSCAN iterations....\n')
    %%
    fprintf(['Adding data to structure....\n'])
    data=struct;
    data.ANCHOR.clusters=ANCHOR_clusters_all;
    data.DBSCAN.clusters = DBSCAN_clusters_all;
    data.CoefRows = CoefRows;
    data.epsilon = epsilon_all;
    data.minpts = minpts_all;
    data.C1 = C1_all;
    data.ANCHOR_coefrowidx= ANCHORidx;
    data.DBSCAN_coefrowidx= DBSCANidx;
    %% REMOVES ORPHANS FROM CACNES ANCHORCLUST
    fprintf(['Starting adoption....\n'])
    Algorithms={'ANCHOR','DBSCAN'};
    A=numel(Algorithms);
    for a = 1:A
        Algo=data.(Algorithms{a});
        clusters=data.(Algorithms{a}).clusters;
        [clusters_unique,data.(Algorithms{a}).ia,data.(Algorithms{a}).ic]=unique(clusters,'rows');
        N = size(clusters_unique,1);
        consistency=false(N,1);
        parclusters_unique=parallel.pool.Constant(clusters_unique);
        DM = parallel.pool.Constant(distance_matrix);
        fprintf(['measuring consistency....\n'])
        parfor n = 1:N
            consistency(n)= check_consistency_fast(DM.Value,parclusters_unique.Value(n,:));
        end
        if any(consistency)
            clusters_unique_consistent=clusters_unique(consistency,:);
            data.(Algorithms{a}).clusters_unique_consistent =clusters_unique_consistent ;
            data.(Algorithms{a}).consistency = consistency;
            U = size(clusters_unique_consistent,1);
            clusters_unique_consistent = parallel.pool.Constant(clusters_unique_consistent);
            clusters_unique_consistent_adopted = zeros(U,S);
            [N,W,C] = deal(zeros(U,1));
            dm = parallel.pool.Constant(distance_matrix);
            fprintf(['adopting lonely isolates....\n'])
            parfor u=1:U
                cl= place_orphans(dm.Value,clusters_unique_consistent.Value(u,:),75);
                clusters_unique_consistent_adopted(u,:)=cl;
                GoodCluster = cl>0;
                b = (cl==cl')&(GoodCluster&GoodCluster');
                W(u) = max(max(dm.Value(b)));
                N(u)=sum(cl>0);
                C(u)=max(cl);
            end
            data.(Algorithms{a}).(fnAlgo{f}).clusters_unique_consistent_adopted=clusters_unique_consistent_adopted;
        end
    end
    fprintf(['Adding data to main structure....\n'])
    SepidermidisATCC12228(d).data=data;
    fprintf(['done with iteration....\n'])
end
%%
fprintf(['saving all data...\n'])
save(['SepidermidisATCC12228.mat'],'SepidermidisATCC12228','-v7.3')