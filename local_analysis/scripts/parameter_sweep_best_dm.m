%% Initialize parallel pool on VM
% distcomp.feature( 'LocalUseMpiexec', false );
% parpool(128,'IdleTimeout',Inf);
% dbstop if error

%% load distance matrices
load data/distance_matrices.mat CacnesDM SepiDM

%% Initialize arrays of clustering parameters
% clustering coefficients
Ucoeficients= .5:.01:1; % cutoffs
Umargins = [1:1:100 200:100:1000]; % distance at which to measure coefficient

% ANCHOR specific
Uneighborhoodr=1:100;

% DBSCAN specific
uminpts=1:3;
uepsilon = 1:1:100;
minnum=1;

%% Index all combinations of parameters for coefficient cutoffs
midx = 1:numel(Umargins);
M = numel(midx);
coef_combinations = [numel(Ucoeficients) numel(Umargins)];
C = prod(coef_combinations);% total possible number of combinations
[a,b]=ind2sub(coef_combinations,1:C);% two dimentional indices
co = Ucoeficients(a)'; % all coefficients
marid = midx(b)'; % all margins

%% Index all combinations of parameters for DBSCAN

iterations = [numel(uepsilon) numel(uminpts)];
Udb = prod(iterations);
[e,m]=ind2sub(iterations,1:Udb);
epsilon =uepsilon(e);
minpts =uminpts(m);

%% Run all instances for both species
% initialize for parallel computing
um = parallel.pool.Constant(Umargins);
for species = 1:2
    % simple logical to switch species without creating a sub-function
    if species ==1
        distance_matrix=SepiDM;
    else
        distance_matrix=CacnesDM;
    end

    %% Initialize parameters for parallel computing

    S = (size(distance_matrix,1));
    coefs_all = zeros(M,S);
    dm=parallel.pool.Constant(distance_matrix);

    %% measuring coefs

    fprintf('measuring clustering coefficients....\n')
    parfor acc = 1:M
        coefs_all(acc,:) = ACcoef(dm.Value,um.Value(acc));
    end

    %% Choose seed isolate

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

    %% Some parameter combinations lead to the same exact sets of isolates.
    % only run unique rows for efficiency.

    fprintf('getting unique rows....\n')
    [UrowsMeetsCutoff,ia,ic]=unique(meets_cutoff_all,'rows');
    CoefRows.UrowsMeetsCutoff=UrowsMeetsCutoff;
    CoefRows.ic=ic;
    CoefRows.ia=ia;
    CoefRows.seed=seed;
    fprintf('done with get_coef_rows....\n')

    %% run ANCHOR for each instance
    fprintf('running ANCHOR....\n')

    % initialize outputs
    R = size(CoefRows.UrowsMeetsCutoff,1);
    Uanchor=numel(Uneighborhoodr);
    [ANCHOR_clusters_all] = deal(zeros(R*Uanchor,S,'int16'));
    [ANCHORidx, C1_all]= deal(zeros(R*Uanchor,1));

    fprintf('starting ANCHOR iterations....\n')
    for r = 1:R % for each unique instance
        idx2 = CoefRows.UrowsMeetsCutoff(r,:)==1;
        dm=parallel.pool.Constant(distance_matrix(idx2,idx2));
        seed=parallel.pool.Constant(CoefRows.seed(idx2));
        idx1 = (((r-1)*Uanchor)+1):(r*Uanchor);
        ANCHORidx(idx1)=r;
        C1_all(idx1)=Uneighborhoodr;
        I=numel(idx1);
        [rclusters_all]=deal(zeros(numel(idx1),sum(idx2),'int16'));
        Uneighborhood=parallel.pool.Constant(Uneighborhoodr);
        for i=1:I % for each neighborhood value
            C1=Uneighborhood.Value(i);
            [anchors,clusters]=(AnchorClust(dm.Value,seed.Value,C1));
            rclusters_all(i,:)=clusters;
        end
        ANCHOR_clusters_all(idx1,idx2)=rclusters_all;
        fprintf([char(string((r/R)*100)) ' percent complete....\n'])
    end
    fprintf('finished ANCHOR iterations....\n')

    %% run DBSCAN

    % initialize
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
    data.clustcoef = co;
    data.cutoffcoef = marid;

    %% pull out consistent instances

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
            consistency(n)= check_consistency(DM.Value,parclusters_unique.Value(n,:));
        end
        if any(consistency)
            clusters_unique_consistent=clusters_unique(consistency,:);
            data.(Algorithms{a}).clusters_unique_consistent =clusters_unique_consistent;
            data.(Algorithms{a}).consistency = consistency;
        end
    end
    fprintf(['done with iteration....\n'])
    fprintf(['saving all data...\n'])
    if species ==1
        SepidermidisATCC12228=data;
        save(['data/SepidermidisATCC12228_clustering_instances.mat.mat'],'SepidermidisATCC12228','-v7.3')
    elseif species==2
        Pacnes_C1=data;
        save(['data/Pacnes_C1_clustering_instances.mat.mat'],'Pacnes_C1','-v7.3')
    end

end