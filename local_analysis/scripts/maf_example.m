function maf_example(unfiltered_mutation_vars,clusters,bool_all,bool_kept)
%%

bool_all = find(bool_all);
bool_kept = find(bool_kept);
bool_clustered=bool_kept(clusters>0);

bool_kept=ismember(bool_all,bool_kept);
bool_clustered=ismember(bool_all,bool_clustered);
%%


% get calls
maNT = single(unfiltered_mutation_vars.maNT_unfiltered(:,bool_all));
maf = single(unfiltered_mutation_vars.maf_unfiltered(:,bool_all));
coverage = single(unfiltered_mutation_vars.coverage_unfiltered(:,bool_all));
Quals = single(unfiltered_mutation_vars.Quals_unfiltered(:,bool_all));
%
Calls = maNT;
% apply first set of filters to Calls_filtered
params.min_qual_for_call = 45;
params.min_maf_for_call = .79;
params.min_cov = 3;
Calls( Quals < params.min_qual_for_call | maf < params.min_maf_for_call | coverage < params.min_cov ) = 0;

P = size(Calls,1);


%% get the NT at each position for each lineage, only where there's one NT that isn't N

%%
proportion_bad_calls = mean(Calls==0,1,'omitnan');
H1=histcounts(proportion_bad_calls(bool_kept&bool_clustered),0:.02:1);
H1=H1./sum(H1);
H2=histcounts(proportion_bad_calls(bool_kept&~bool_clustered),0:.02:1);
H2=H2./sum(H2);
% if any(proportion_bad_calls(~bool_kept))
%     H3=histcounts(proportion_bad_calls(~bool_kept),0:.02:1);
%     H3=H3./sum(H3);
%     b=bar(([H1 ; H2; H3])','stacked','FaceColor','flat')
%     b(1).CData=[0 0 1]
%     b(2).CData=[1 0 0]
%     b(3).CData=[.85 .85 .85]
% else
    b=bar(([H1 ; H2])','stacked','FaceColor','flat')
    b(1).CData=[0 0 1]
    b(2).CData=[1 0 0]
% end
bin_right_edges=.02:.02:1;
xticks([0 5:5:50])
xticklabels([0 .1:.1:1])
xlabel('proportion of positions with bad calls')
ylabel('proportion samples')
[~,pv,~] = kstest2(proportion_bad_calls(bool_kept&bool_clustered),proportion_bad_calls(bool_kept&~bool_clustered))
%%
