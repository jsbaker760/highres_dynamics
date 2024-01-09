%% MAIN SCRIPT FOR PARALLEL EVOLUTION ANALYSIS
function annotations_mutgenes_all = evo_parallel_main_wrapper(annotations_all,mutation_spectrum,species)
% Last updated: Arolyn, 2022.05.06 for JSB Acera S epi analysis





%% Nucleotides: 1=A, 2=T, 3=C, 4=G





NTs='ATCG';





%% create directory for figures





dir_summ = 'figs_summary_parallel_evolution';
if ~exist(dir_summ,'dir')
    mkdir(dir_summ)
end





%% Load positions masked in SNP calling





positions_to_mask=[]; % empty for now





%% GENERAL PARAMETERS FOR PARALLEL EVO ANALYSIS





% Find number of genes under parallel evolution detected at these false
% discovery rates (Benhamini-Hochberg FDR):

fdr_to_test = [ 0 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 ];

% Compute dN/dS for this many genes with the lowest Poisson p-values
top_num_genes_to_test = [ 2, 5, 10 ];

% Minimum number of mutations to evaluate set
min_num_muts = 3;





%% ANALYZE MUTATIONS BY LINEAGE





Ulineages = unique([annotations_all(:).cluster]);

lineage_names = arrayfun(@(x) {char(strjoin([species "clade" string(x)],'_'))}, Ulineages);

% Initialize objects for storing results
annotations_mutgenes_all = struct;

for lin_index=1:numel(lineage_names)

    this_lineage = lineage_names{lin_index};

    annotation_full_everything = annotations_all( [annotations_all.cluster] == lin_index );
    data_set_name = this_lineage;
    
    % Where to find the reference genome on the Lieberman Lab Dropbox:
    dir_ref_genome = [ 'data/Assembly/clade_coassemblies/clades_filtered_contigs/'  this_lineage];
    [~, GenomeLength, ~, ~] = genomestats(dir_ref_genome);

    if length( annotation_full_everything ) >= min_num_muts
        % Call parallel evolution function
        [ annotations_mutgenes ] ...
            = evo_parallel_main( data_set_name, annotation_full_everything, ...
            dir_ref_genome, GenomeLength,positions_to_mask, mutation_spectrum.mut_spec_prob, ...
            fdr_to_test, top_num_genes_to_test );
        % Append to annotations_mutgenes_all
        if numel(fieldnames(annotations_mutgenes_all))==0
            annotations_mutgenes_all = annotations_mutgenes;
        else
            annotations_mutgenes_all = [ annotations_mutgenes_all, annotations_mutgenes ];
        end
    else
        fprintf(1,['Insufficient mutations in Lineage ' this_lineage ' (n_muts=' num2str(length(annotation_full_everything)) ')' '\n'])
    end
    
end

