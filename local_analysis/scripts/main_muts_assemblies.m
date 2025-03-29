function [annotations_geneclusters,all_proteins_info_assemblies,MutationSpectrum_obs,clade_expected_pnonsyn_mean]=main_muts_assemblies(annotations_all,HomologousClustersAssemblies,clade_ref_dir)


NTs='ATCG';

%% Load data about all within-lineage mutations

% Number of mutations
num_muts_all = numel(annotations_all);

%% Compute observed mutation spectrum


% mutationmatrix ATCG -> ATCG matrix; gathered into mut_observed (6 entry vector)
% mut_observed tally of mutations of 6 types
% prob_nonsyn is the calculated probability of a mutation being N

% Compute mutation spectrum based on observed SNPs
[ mutationmatrix, mut_observed, typecounts, prob_nonsyn ] = ...
    mutation_spectrum_module( annotations_all, NTs);
mut_spec_prob = mut_observed/sum(mut_observed);
mut_spec_names = { 'AT/TA', 'AC/TG', 'AG/TC', 'GC/CG', 'GT/CA', 'GA/CT' };

% from mutation_spectrum_module via div_matrix2_6types
% AT, TA % transversion
% AC, TG % transversion
% AG, TC % transition
% GC, CG % transversion
% GT, CA % transversion
% GA, CT % transition

MutationSpectrum_obs=struct;
MutationSpectrum_obs.mutationmatrix=mutationmatrix;
MutationSpectrum_obs.mut_observed=mut_observed;
MutationSpectrum_obs.typecounts=typecounts;
MutationSpectrum_obs.prob_nonsyn=prob_nonsyn;
MutationSpectrum_obs.mut_spec_names=mut_spec_names;
MutationSpectrum_obs.mut_spec_prob=mut_spec_prob;

%% Load information about genomes

% the lineages which have mutation

ulineages=unique([annotations_all.cluster]);
% Genome lengths, expected p_nonsyn
num_clades = numel(ulineages);
clade_genome_lengths = zeros(num_clades,1);
clade_expected_pnonsyn = zeros(num_clades,1);
parfor c=1:num_clades
    clade_number = ulineages(c);
    fprintf(1,['Clade ' num2str(clade_number) '...\n'])
    % Genome length
    dir_ref_genome = [ clade_ref_dir num2str(clade_number) ];
    [~, GenomeLength, ~, ~] = genomestats(dir_ref_genome);
    clade_genome_lengths(c) = GenomeLength;
    % Expected p_nonsyn
    positions_to_mask = [];
    cdsindex_proteins = get_coding_genenums( dir_ref_genome, positions_to_mask ); % avoid genes in regions masked in SNP calling
    clade_expected_pnonsyn(c) = compute_expected_dnds( dir_ref_genome, mut_spec_prob, cdsindex_proteins );
end
% the expected probability of a non-syn mutation
clade_expected_pnonsyn_mean = mean(clade_expected_pnonsyn);

%% show that p(nonsyn) is very similar across genomes

figure(10)
clf(10)
hold on
box on
histogram( clade_expected_pnonsyn,BinWidth=.001 )
xlabel('probability random mutation is nonsyn')
ylabel('number of clades')
set(gca,'FontSize',14)
hold off


%% Number of mutations per clade (0 where there are no mutations)
num_muts_per_clade = zeros(1,max(ulineages));
for c=1:num_clades
    num_muts_per_clade(c) = sum( c==[annotations_all.cluster] );
end


%% get relavent clade numbers from Roary data

ulineages=unique(HomologousClustersAssemblies.CladeNumber);

%% get all_proteins_info, a structure with information for every protein in each genome assembly (one for each clade)

all_proteins_info_assemblies = get_all_protein_info(ulineages,MutationSpectrum_obs.mut_spec_prob,clade_ref_dir);

%% get the locustags from each protein for each assembly

all_proteins_info_assemblies=struct2table(all_proteins_info_assemblies);
locustag_all_proteins_info=string(all_proteins_info_assemblies.locustag);

%% add homologous cluster to all protein info

locustag_homologs= string([HomologousClustersAssemblies.GeneTag]);
homolog_numbers = ([HomologousClustersAssemblies.HomologueNumber]);

all_proteins_info_gene_cluster = arrayfun(@(x) max([0 homolog_numbers(locustag_homologs==x)]), locustag_all_proteins_info);
all_proteins_info_assemblies.gene_cluster=all_proteins_info_gene_cluster;
all_proteins_info_assemblies=table2struct(all_proteins_info_assemblies);

%% dd homologous cluster to annotations_all

table_annotations_all = struct2table(annotations_all);
locustags = table_annotations_all.locustag;
locustags(cellfun(@isempty,locustags))={''}; % mutations which are not in genes
locustags=string(locustags);
mutation_homologs = arrayfun(@(x) max([0 HomologousClustersAssemblies.HomologueNumber(HomologousClustersAssemblies.GeneTag==x)]), locustags);

% turn it back into a structure
table_annotations_all.HomologousCluster = mutation_homologs;
annotations_all=table2struct(table_annotations_all);


%% Observed mutations aggregated into each homologous protein cluster

annotations_geneclusters = struct;
num_geneclusters_with_muts = 0;
unique_mutated_geneclusters = unique(mutation_homologs(mutation_homologs>0)); % unique mutated homologous clusters
num_gene_clusters=numel(unique_mutated_geneclusters);

% iterate through each

for i=1:num_gene_clusters % for each unique homolgous cluster

    next_genecluster=unique_mutated_geneclusters(i);
    genecluster_mut_idx = find(mutation_homologs==next_genecluster);
    num_geneclusters_with_muts = num_geneclusters_with_muts + 1;

    % Observed number of mutations in this gene cluster across all clades
    next_num_muts_obs = numel(genecluster_mut_idx);
    fprintf( [ 'Gene cluster: ' num2str(next_genecluster) '. Num muts: ' num2str(next_num_muts_obs) '\n' ] )

    % Expected number of mutations in this gene cluster across all claders
    geneclusters_protein_idx = find( next_genecluster==[all_proteins_info_assemblies.gene_cluster] ); % index of protein across all proteins in all assemblies
    next_protein_clades = [all_proteins_info_assemblies(geneclusters_protein_idx).clade]; % the lineages of these proteins
    next_protein_probs =  [all_proteins_info_assemblies(geneclusters_protein_idx).prob_gene_mutated]; % the probabilities they are mutated
    next_num_muts_expected = sum( num_muts_per_clade(next_protein_clades).*next_protein_probs );

    % Poisson pvalue for the observed number of mutations
    next_pval_poisson = poisscdf( next_num_muts_obs-1, next_num_muts_expected, 'upper' ); % poisscdf( x, labmda, 'upper' ) computes pval for >x not >=x

    % Observed N/S
    mutation_types = [annotations_all(genecluster_mut_idx).type];
    obs_num_N = sum(mutation_types=='N');
    obs_num_S = sum(mutation_types=='S');
    [ p_nonsyn_observed, CI_nonsyn ] = binofit( obs_num_N, obs_num_N + obs_num_S, 0.05 );
    NS_observed = p_nonsyn_observed/(1-p_nonsyn_observed);
    NS_obs_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) );
    NS_obs_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) );

    % Expected N/S
    next_genomes = [annotations_all(genecluster_mut_idx).cluster];
    next_genenums = [annotations_all(genecluster_mut_idx).gene_num];
    p_nonsyn_expected_list = zeros(1,numel(next_genomes));

    for j=1:numel(next_genomes)
        p_nonsyn_expected_list(j) = compute_expected_dnds( [ clade_ref_dir num2str(next_genomes(j)) ], mut_spec_prob, next_genenums(j) );
    end
    NS_expected_list = p_nonsyn_expected_list./(1-p_nonsyn_expected_list);
    NS_expected = mean(NS_expected_list);

    % Observed dN/dS
    dNdS_observed = NS_observed/NS_expected;
    dNdS_obs_lower = NS_obs_lower/NS_expected;
    dNdS_obs_upper = NS_obs_upper/NS_expected;

    % Annotation and protein sequence - just grab the first one
    if numel(genecluster_mut_idx)>1
        next_annotation = annotations_all(genecluster_mut_idx(1)).protein;
        next_nt_sequence = annotations_all(genecluster_mut_idx(1)).Sequence;
        next_aa_sequence = annotations_all(genecluster_mut_idx(1)).translation;
        next_gene_length = abs(annotations_all(genecluster_mut_idx(1)).loc2-annotations_all(genecluster_mut_idx(1)).loc1)+1;
    else
        next_annotation = annotations_all(genecluster_mut_idx).protein;
        next_nt_sequence = annotations_all(genecluster_mut_idx(1)).Sequence;
        next_aa_sequence = annotations_all(genecluster_mut_idx).translation;
        next_gene_length = abs(annotations_all(genecluster_mut_idx).loc2-annotations_all(genecluster_mut_idx).loc1)+1;
    end
    %%
    % Record findings
    annotations_geneclusters(num_geneclusters_with_muts).genecluster = next_genecluster;
    annotations_geneclusters(num_geneclusters_with_muts).num_muts_obs = next_num_muts_obs;
    annotations_geneclusters(num_geneclusters_with_muts).num_muts_expected = next_num_muts_expected;
    annotations_geneclusters(num_geneclusters_with_muts).pval_poisson = next_pval_poisson;
    annotations_geneclusters(num_geneclusters_with_muts).num_N = obs_num_N;
    annotations_geneclusters(num_geneclusters_with_muts).num_S = obs_num_S;
    annotations_geneclusters(num_geneclusters_with_muts).NS_obs = NS_observed;
    annotations_geneclusters(num_geneclusters_with_muts).NS_obs_lower = NS_obs_lower;
    annotations_geneclusters(num_geneclusters_with_muts).NS_obs_upper = NS_obs_upper;
    annotations_geneclusters(num_geneclusters_with_muts).NS_expected = NS_expected;
    annotations_geneclusters(num_geneclusters_with_muts).dNdS_obs = dNdS_observed;
    annotations_geneclusters(num_geneclusters_with_muts).dNdS_obs_lower = dNdS_obs_lower;
    annotations_geneclusters(num_geneclusters_with_muts).dNdS_obs_upper = dNdS_obs_upper;
    annotations_geneclusters(num_geneclusters_with_muts).protein = next_annotation;
    annotations_geneclusters(num_geneclusters_with_muts).sequence_ex = next_nt_sequence;
    annotations_geneclusters(num_geneclusters_with_muts).translation_ex = next_aa_sequence;
    annotations_geneclusters(num_geneclusters_with_muts).gene_length_ex = next_gene_length;
end


%% end of function
end