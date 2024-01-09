function [annotations_geneclusters,assemblies_all_proteins_info]=parevo_homologues(annotations_all,MutationSpectrum,Roary,species_str)
    %% Set up directories and environment
    
clade_ref_dir=['data/Assembly/clade_coassemblies/clades_filtered_contigs/' species_str '_clade_'];

dir_summ = 'figs_summary';
if ~exist(dir_summ,'dir')
    mkdir(dir_summ)
end



%% get relavent clade numbers from Roary data 

ulineages=unique(Roary.CladeNumber);

%% get all_proteins_info, a structure with information for every protein in each genome assembly (one for each clade)

assemblies_all_proteins_info = get_all_protein_info(ulineages,MutationSpectrum.mut_spec_prob,clade_ref_dir);
% parse the locustags from all_protein_info
assemblies_all_locus_tags = {assemblies_all_proteins_info(:).locustag};
assemblies_all_locus_tags(cellfun(@isempty,assemblies_all_locus_tags))={''};
assemblies_all_locus_tags=string(assemblies_all_locus_tags);

%% Add homologues cluster to info for each protein in each assembly


assemblies_all_homologue_numbers = arrayfun(@(x) {Roary.HomologueNumber(x==Roary.GeneTag)} , assemblies_all_locus_tags);
assemblies_all_homologue_numbers(cellfun(@isempty,assemblies_all_homologue_numbers))={0};
% put the homologue cluster for each protein
[assemblies_all_proteins_info.HomologousCluster]=assemblies_all_homologue_numbers{:};

%% NTs
% Nucleotides: 1=A, 2=T, 3=C, 4=G

NTs='ATCG';

%% number of mutations in annotations_all

% Number of mutations
N = numel(annotations_all);
L=numel(ulineages);


%% get p(nonsynonymous) and length for each genome
% Initialize genome info for lineages
% Genome lengths, expected p_nonsyn
% Number of mutations per clade
uclade_genome_lengths = zeros(L,1);
uclade_expected_pnonsyn = zeros(L,1);
uclades_num_muts = zeros(1,L);

for l=1:L
    fprintf(1,['Clade ' num2str(ulineages(l)) '...\n'])
    % Genome length
    dir_ref_genome = [clade_ref_dir char(string(ulineages(l)))];
    [~, GenomeLength, ~, ~] = genomestats(dir_ref_genome);
    uclade_genome_lengths(l) = GenomeLength;
    % Expected p_nonsyn
    cdsindex_proteins = get_coding_genenums( dir_ref_genome); % avoid genes in regions masked in SNP calling
    uclade_expected_pnonsyn(l) = compute_expected_dnds( dir_ref_genome, MutationSpectrum.mut_spec_prob, cdsindex_proteins );
    uclades_num_muts(l) = sum( ulineages(l)==[annotations_all.cluster] );
end



%% Get the homologous cluster for each annotation

anno_locustags={annotations_all(:).locustag} ;
anno_locustags(cellfun(@isempty,anno_locustags))={''};
anno_locustags=string(anno_locustags);
% add it where these is a locustag

%%
assemblies_all_homologue_numbers=vertcat([assemblies_all_homologue_numbers{:}]);
%%
annotation_homologous_clusters=arrayfun(@(x) max([0 assemblies_all_homologue_numbers(x==assemblies_all_locus_tags)]), anno_locustags);

%% Observed mutations aggreated into each protein cluster

UniqueMutatedHomologues= unique(annotation_homologous_clusters(annotation_homologous_clusters>0));
annotations_geneclusters = struct;
for g=1:numel(UniqueMutatedHomologues)
    ThisProteinCluster=UniqueMutatedHomologues(g);
    idxAnnotationThisHomologue = find(annotation_homologous_clusters==ThisProteinCluster);
    if ~isempty(idxAnnotationThisHomologue)
        % counting the number of geneclusters with mutations

        % Observed number of mutations in this gene cluster across all clades
        num_muts_obs_gene = numel(idxAnnotationThisHomologue);
        fprintf( [ 'Gene cluster: ' num2str(ThisProteinCluster) '. Num muts: ' num2str(num_muts_obs_gene) '\n' ] )
        % Expected number of mutations in this gene cluster across all claders
        idx_homologue_in_assemblies = find( ThisProteinCluster==assemblies_all_homologue_numbers);
        names_clades_with_homologue = [assemblies_all_proteins_info(idx_homologue_in_assemblies).clade];
        protein_probs_clades_with_homologue =  [assemblies_all_proteins_info(idx_homologue_in_assemblies).prob_gene_mutated];
        % index in the size of ulineages
        idx_lineages=arrayfun(@(x) find(x==ulineages), names_clades_with_homologue);
        protein_probs_clades_with_homologue=protein_probs_clades_with_homologue(idx_lineages);
        num_muts_expected = sum( uclades_num_muts(idx_lineages).*protein_probs_clades_with_homologue);
        
      
        % Poisson pvalue for the observed number of mutations
        next_pval_poisson = poisscdf( num_muts_obs_gene-1, num_muts_expected, 'upper' ); % poisscdf( x, labmda, 'upper' ) computes pval for >x not >=x

        % Observed N/S
        next_types = [annotations_all(idxAnnotationThisHomologue).type];
        next_num_N = sum(next_types=='N');
        next_num_S = sum(next_types=='S');
        [ p_nonsyn_observed, CI_nonsyn ] = binofit( next_num_N, next_num_N + next_num_S, 0.05 );
        NS_observed = p_nonsyn_observed/(1-p_nonsyn_observed);
        NS_obs_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) );
        NS_obs_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) );

        % Expected N/S
        next_genomes = [annotations_all(idxAnnotationThisHomologue).cluster];
        next_genenums = [annotations_all(idxAnnotationThisHomologue).gene_num];
        p_nonsyn_expected_list = zeros(1,numel(next_genomes));

        for i=1:numel(next_genomes)
            % #TODO: rewrite so genomes don't have to be loaded live
            p_nonsyn_expected_list(i) = compute_expected_dnds( [ clade_ref_dir num2str(next_genomes(i)) ], MutationSpectrum.mut_spec_prob, next_genenums(i) );
        end

        NS_expected_list = p_nonsyn_expected_list./(1-p_nonsyn_expected_list);
        NS_expected = mean(NS_expected_list);


        % Observed dN/dS
        dNdS_observed = NS_observed/NS_expected;
        dNdS_obs_lower = NS_obs_lower/NS_expected;
        dNdS_obs_upper = NS_obs_upper/NS_expected;

        % Annotation and protein sequence - just grab the first one
        if numel(idxAnnotationThisHomologue)>1
            next_annotation = annotations_all(idxAnnotationThisHomologue(1)).protein;
            next_nt_sequence = annotations_all(idxAnnotationThisHomologue(1)).Sequence;
            next_aa_sequence = annotations_all(idxAnnotationThisHomologue(1)).translation;
            next_gene_length = abs(annotations_all(idxAnnotationThisHomologue(1)).loc2-annotations_all(idxAnnotationThisHomologue(1)).loc1)+1;
        else
            next_annotation = annotations_all(idxAnnotationThisHomologue).protein;
            next_nt_sequence = annotations_all(idxAnnotationThisHomologue(1)).Sequence;
            next_aa_sequence = annotations_all(idxAnnotationThisHomologue).translation;
            next_gene_length = abs(annotations_all(idxAnnotationThisHomologue).loc2-annotations_all(idxAnnotationThisHomologue).loc1)+1;
        end



        %%
        % Record findings
        annotations_geneclusters(g).genecluster = ThisProteinCluster;
        annotations_geneclusters(g).num_muts_obs = num_muts_obs_gene;
        annotations_geneclusters(g).num_muts_expected = num_muts_expected;
        annotations_geneclusters(g).pval_poisson = next_pval_poisson;
        annotations_geneclusters(g).num_N = next_num_N;
        annotations_geneclusters(g).num_S = next_num_S;
        annotations_geneclusters(g).NS_obs = NS_observed;
        annotations_geneclusters(g).NS_obs_lower = NS_obs_lower;
        annotations_geneclusters(g).NS_obs_upper = NS_obs_upper;
        annotations_geneclusters(g).NS_expected = NS_expected;
        annotations_geneclusters(g).dNdS_obs = dNdS_observed;
        annotations_geneclusters(g).dNdS_obs_lower = dNdS_obs_lower;
        annotations_geneclusters(g).dNdS_obs_upper = dNdS_obs_upper;
        annotations_geneclusters(g).protein = next_annotation;
        annotations_geneclusters(g).sequence_ex = next_nt_sequence;
        annotations_geneclusters(g).translation_ex = next_aa_sequence;
        annotations_geneclusters(g).gene_length_ex = next_gene_length;
        annotations_geneclusters(g).annotations_indices= idxAnnotationThisHomologue;
    end
end



