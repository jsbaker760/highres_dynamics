
CacnesClusters = unique([Pacnes_C1.ANCHOR.clusters_unique_consistent ; Pacnes_C1.DBSCAN.clusters_unique_consistent],'rows');
CacnesClusters((sum(CacnesClusters>0,2)<1800)|max(CacnesClusters,[],2)<96|max(CacnesClusters,[],2)>150,:)=[];


%%

SepiClusters = unique([SepidermidisATCC12228.ANCHOR.clusters_unique_consistent ; SepidermidisATCC12228.DBSCAN.clusters_unique_consistent],'rows');
SepiClusters(sum(SepiClusters>0,2)<1600|max(SepiClusters,[],2)<100|max(SepiClusters,[],2)>200,:)=[];

%%
CacnesClustersNew = fix_it(CacnesDM,CacnesClusters);
save('CacnesClustersNew.mat','CacnesClustersNew','-v7.3')
%%
% load('GoodDMs.mat', 'SepiDM')
SepiClustersNew = fix_it(SepiDM,SepiClusters);
%%
save('SepiClustersNew.mat','SepiClustersNew','-v7.3')
%%


function newclusters = fix_it(DM,clusters_old)

[U,S] = size(clusters_old);
newclusters = zeros(U,S);
clusters_old=parallel.pool.Constant(clusters_old);
DM = parallel.pool.Constant(DM);

showTimeToCompletion; startTime=tic;
p = parfor_progress(U);


parfor u=1:U

    n= place_orphans(DM.Value,clusters_old.Value(u,:),75);
    newclusters(u,:)=n;
    p = parfor_progress;
    showTimeToCompletion( p/100, [], [], startTime );
end

end