function [f1,f2,f3]=plot_evolution_homologous_genes(annotations_geneclusters,probnonsyn)

%% get relavent variables from structure for Benjamini-Hochberg analysis

InfSubstitute=200;

%% for sorting genes by FDR
pval_poisson_list = [annotations_geneclusters.pval_poisson];

% for dN/dS on genes by Nmutations/Mutation density
obs_num_muts_per_gene=[annotations_geneclusters.num_muts_obs];
gene_lengths=[annotations_geneclusters.gene_length_ex];
obs_num_mutsper1kb_per_gene=(obs_num_muts_per_gene./gene_lengths)*1000;

%% list of proportions (percent/100) to compare in FDR

fdr_bh_list = [.01 .05 .1 .5 1 10 100];

%% plot dN/dS versus BHFDR

f1=figure;hold on;

for i=1:numel(fdr_bh_list)
    fdr_bh = fdr_bh_list(i); % false discovery rate
    psort = sort(pval_poisson_list,'ascend');
    alpha = ((1:numel(psort))/numel(psort))*fdr_bh;

    crit=max(psort(psort<=alpha));
    if ~isempty(crit)
        genes_of_interest_bool = ( pval_poisson_list<=crit );
        num_N = sum([annotations_geneclusters(genes_of_interest_bool).num_N]);
        num_S = sum([annotations_geneclusters(genes_of_interest_bool).num_S]);

        [dnds_set,dnds_set_lower,dnds_set_upper]=dnds_homologues(num_N,num_S,probnonsyn);

        bar(i,log2(dnds_set),'BaseValue',0,'FaceColor',[.85 .85 .85],'LineWidth',1)
        plot([i i], [log2(dnds_set_lower) log2(dnds_set_upper)],'Color','black','LineWidth',1)
    end
end

plot([0 numel(fdr_bh_list)+1],[0 0],'Color','black','LineWidth',1)

xticks([1:numel(fdr_bh_list)])
xticklabels(fdr_bh_list)
xlabel('Benjamini-Hochberg FDR')
xlim([0 numel(fdr_bh_list)+1])
h=gca;
h.LineWidth=1;

yt=-1:1:1;
realvals=2.^yt
yticks(yt);
yticklabels(realvals)
ylim([-1 1]);
ylabel('dN/dS')

%% dN/dS for genes above a certain number of mutation

% Filters
list_nmuts_thresholds = [1 2 5 10];
L=numel(list_nmuts_thresholds);

% initialize
[dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset]=deal(zeros(L,1));

f2=figure;hold on;

for i=1:L
    % true where gene is mutated at high enough density and has at least 2
    % mutations
    genes_of_interest = obs_num_muts_per_gene >= list_nmuts_thresholds(i);
    %
    num_N = sum([annotations_geneclusters(genes_of_interest).num_N]);
    num_S = sum([annotations_geneclusters(genes_of_interest).num_S]);

    [dNdS_by_geneset(i), dNdS_lower_by_geneset(i), dNdS_upper_by_geneset(i)]=dnds_homologues(num_N,num_S,probnonsyn);
    
    bar(i,log2(dNdS_by_geneset(i)),'BaseValue',0,'FaceColor',[.85 .85 .85],'LineWidth',1)
    plot([i i], [log2(dNdS_lower_by_geneset(i)) log2(dNdS_upper_by_geneset(i))],'Color','black','LineWidth',1)

end

plot([0 L+1],[0 0],'Color','black','LineWidth',1)


xticks([1:numel(list_nmuts_thresholds)])
xticklabels(list_nmuts_thresholds)
xlabel('Mutation threshhold (nmuts/gene)')

yt=-1:1:1;
realvals=2.^yt
yticks(yt);
yticklabels(realvals)
ylabel('dN/dS')
ylim([-1 1])
xlim([0 numel(list_nmuts_thresholds)+1])
h=gca;
h.LineWidth=1;

%% dN/dS for genes above a certain mutational density

% Filters
min_num_muts_per_gene = 2; % only do density for multiply mutated genes
list_density_thresholds = [ 0 1 2 5 ];
L=numel(list_density_thresholds);

% initialize
[dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset]=deal(zeros(L,1));

f3=figure;hold on;

for i=1:L
    % true where gene is mutated at high enough density and has at least 2
    % mutations
    if list_density_thresholds(i)==0
        genes_of_interest=obs_num_muts_per_gene>0; % turn it off to have comparison of all mutated genes
    else
        genes_of_interest = obs_num_muts_per_gene >= min_num_muts_per_gene & obs_num_mutsper1kb_per_gene >= list_density_thresholds(i);
    end
    %
    num_N = sum([annotations_geneclusters(genes_of_interest).num_N]);
    num_S = sum([annotations_geneclusters(genes_of_interest).num_S]);

    [dNdS_by_geneset(i), dNdS_lower_by_geneset(i), dNdS_upper_by_geneset(i)]=dnds_homologues(num_N,num_S,probnonsyn);
    
    bar(i,log2(dNdS_by_geneset(i)),'BaseValue',0,'FaceColor',[.85 .85 .85],'LineWidth',1)
    plot([i i], [log2(dNdS_lower_by_geneset(i)) log2(dNdS_upper_by_geneset(i))],'Color','black','LineWidth',1)

end

plot([0 L+1],[0 0],'Color','black','LineWidth',1)


xticks([1:numel(list_density_thresholds)])
xl=string(list_density_thresholds);
xl(1)="all mutated genes"
xticklabels(xl)
xlabel('Mutation density threshhold (muts/kb)')

yt=-1:1:1;
realvals=2.^yt
yticks(yt);
yticklabels(realvals)
ylabel('dN/dS')
ylim([-1 1])
xlim([0 numel(list_density_thresholds)+1])
h=gca;
h.LineWidth=1;


% end of function
end
%% compute dN/dS for a subset of homologues


function [dNdS, dNdS_lower, dNdS_upper]=dnds_homologues(num_muts_N,num_muts_S,probnonsyn_expected)

%% observed

% Compute probability of a nonsynonymous mutation with 95% CIs using binofit
% [ p, CI ] = binofit( x, n, alpha )
% p = max likelihood prob of succcess
% CI = confidence interval (100(1-alpha))
% x = number of successes
% n = number of trials
% alpha = to specify confidence interval (optional)

% observed ratio of N to S mutations
ns_observed=num_muts_N/num_muts_S;

% confidence interval of observation
x = num_muts_N;
n = num_muts_N + num_muts_S;
alpha = 0.05;
[ ~, CI_nonsyn ] = binofit( x, n, alpha );

% expected N/S from a neutral model
ns_expected_local = probnonsyn_expected/(1-probnonsyn_expected);

% Compute dN/dS
dNdS =  ns_observed/ns_expected_local; % relative to neutral model
dNdS_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / ns_expected_local;
dNdS_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / ns_expected_local;



end