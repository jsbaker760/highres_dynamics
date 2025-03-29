function plot_refpositions_shared(unfiltered_mutation_vars,clusters,bool)

% get calls 
maNT = single(unfiltered_mutation_vars.maNT_unfiltered(:,bool));
maf = single(unfiltered_mutation_vars.maf_unfiltered(:,bool));
coverage = single(unfiltered_mutation_vars.coverage_unfiltered(:,bool));
Quals = single(unfiltered_mutation_vars.Quals_unfiltered(:,bool));
%
Calls = maNT;
% apply first set of filters to Calls_filtered
params.min_qual_for_call = 40;
params.min_maf_for_call = .79;
params.min_cov = 3;
Calls( Quals < params.min_qual_for_call | maf < params.min_maf_for_call | coverage < params.min_cov ) = 0;

P = size(Calls,1);
%% get the NT at each position for each lineage, only where there's one NT that isn't N

uclusters = unique(clusters(clusters>0));
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


%%
has_nt = lineages_calls>0;
%%
T = 100;
NSharedPos = zeros(T,U);
Ncompared = repmat(1:U,T,1);
%
[~,perms_all] = sort(rand(T*U,U),2);
for u=1:U
    perms=perms_all(((T*(u-1))+1):((T*u)),1:u);
    l_all = has_nt(:,perms(:));
    l_all = reshape(l_all,[],T,u);
    universal=sum(~any(~l_all,3),1);
    NSharedPos(:,u)=universal;
    u/U
end
%

scatter(Ncompared,NSharedPos,10,'filled','MarkerFaceColor','black')
