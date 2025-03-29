function snp_distance_distortion(unfiltered_mutation_vars,clusters,bool)
% reduce clusters to only clustered isolates
clustered = clusters>0;
bool = find(bool);
bool = bool(clustered);
clusters=clusters(clustered);


% get calls 
maNT = single(unfiltered_mutation_vars.maNT_unfiltered(:,bool));
maf = single(unfiltered_mutation_vars.maf_unfiltered(:,bool));
coverage = single(unfiltered_mutation_vars.coverage_unfiltered(:,bool));
Quals = single(unfiltered_mutation_vars.Quals_unfiltered(:,bool));
Calls = maNT;

% apply first set of filters to Calls_filtered
params.min_qual_for_call = 45;
params.min_maf_for_call = .79;
params.min_cov = 3;
Calls( Quals < params.min_qual_for_call | maf < params.min_maf_for_call | coverage < params.min_cov ) = 0;

P = size(Calls,1);


% get the NT at each position for each lineage, only where there's one NT that isn't N

uclusters = unique(clusters);
U = numel(uclusters);
lineages_calls = zeros(P,U);

for u = 1:U
    calls_lineage = zeros(P,1);
    calls_isolates = Calls(:,clusters==uclusters(u));
    mode_nt = single(calls_isolates);
    mode_nt(mode_nt==0)=NaN;
    mode_nt=mode(mode_nt,2);
    mode_nt(isnan(mode_nt))=0;
    n_different_alleles = sum(diff(sort(calls_isolates,2),1,2)~=0,2)+1;
    any_n = sum(calls_isolates==0,2)>0;
    only_one_nt = n_different_alleles==1|(n_different_alleles==2&any_n);
    calls_lineage(only_one_nt)=mode_nt(only_one_nt);
    lineages_calls(:,u)=calls_lineage;
end



% all pos

Calls=lineages_calls;
numsamples=size(Calls,2);
goodcall = (Calls>0);

distance_matrix_allpos=(single(zeros(numsamples)));
tic;
for s=1:numsamples
    sample_calls = repmat(Calls(:,s),1,numsamples);
    distance_matrix_allpos(s,:)=sum(Calls~=sample_calls & goodcall & sample_calls>0);
end
toc

%% core pos

proportion_call = mean(Calls>0,2);

Calls=lineages_calls(proportion_call==1,:);
numsamples=size(Calls,2);
goodcall = (Calls>0);

distance_matrix_corepos=(single(zeros(numsamples)));
tic;
for s=1:numsamples
    sample_calls = repmat(Calls(:,s),1,numsamples);
    distance_matrix_corepos(s,:)=sum(Calls~=sample_calls & goodcall & sample_calls>0);
end
toc

%% plot
same_cluster = uclusters==uclusters';
different_cluster = uclusters~=uclusters';

subplot(2,1,1)
scatter(distance_matrix_corepos(different_cluster),distance_matrix_allpos(different_cluster)./distance_matrix_corepos(different_cluster),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','black')
hold on
plot([0 10^5],[1 1],'LineWidth',1,'LineStyle','--','Color','red')
xlim([0 50000])
ylim([0 25])
xlabel('core genome SNP positions')
ylabel('fold-increase in SNPs')
pbaspect([1 1 1])
subplot(2,1,2)
scatter(distance_matrix_corepos(different_cluster),distance_matrix_allpos(different_cluster)./distance_matrix_corepos(different_cluster),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','black')
xlim([0 1000])
ylim([0 5])
xlabel('core genome SNP positions')
ylabel('fold-increase in SNPs')
pbaspect([1 1 1])