function [ annotations_mutgenes ] = evo_parallel_main( data_set_name, annotation_full_everything,dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, fdr_to_test, top_num_genes_to_test )





fprintf(1,['Analyzing data set: ' data_set_name '...\n'])





%% DIRECTORY SETUP
fprintf(1,['Setting up directories...' '\n'])





% Where to save results
dir_results = ['figs_parallelevo_' data_set_name];
if ~exist( dir_results, 'dir')
    mkdir( dir_results )
end

dir_geneinfo = [ dir_results '/' 'gene_info' ];
if ~exist( dir_geneinfo, 'dir')
    mkdir( dir_geneinfo )
end





%% I. OBSERVED and EXPECTED NUMBER OF GENES WITH MULTIPLE MUTATIONS

%% Get gene numbers for proteins only (exclude tRNA and rRNA)
fprintf(1,['Getting genome information...' '\n'])





cdsindex_proteins = get_coding_genenums( dir_ref_genome ); % avoid genes in regions masked in SNP calling
num_proteins_in_genome = numel( cdsindex_proteins );

% looks good
p_nonsyn_expected_allproteins = compute_expected_dnds( dir_ref_genome, mut_spec_prob, cdsindex_proteins );
probstop_expected_allproteins = compute_probstop( dir_ref_genome, mut_spec_prob, cdsindex_proteins );
NS_expected_allproteins = p_nonsyn_expected_allproteins/(1-p_nonsyn_expected_allproteins);
STOP_expected_allproteins = probstop_expected_allproteins/(1-probstop_expected_allproteins);





%% Observed genes with multiple mutations
fprintf(1,['Finding mutated genes...' '\n'])





[ obs_num_muts_per_gene, obs_mutated_genes_genenums, obs_genes_lengths ] = ...
    compute_observed_parallelevo( annotation_full_everything, cdsindex_proteins ); % sim_gene_genenums used to exclude rRNA, tRNA
obs_num_mutsper1kb_per_gene = 1000*obs_num_muts_per_gene./obs_genes_lengths;
num_genes_mutated = numel( obs_mutated_genes_genenums );
obs_num_muts_total = sum( obs_num_muts_per_gene );





%% Expected genes with multiple mutations (simulate neutral model)
fprintf(1,['Simulating mutated genes...' '\n'])





num_muts_to_sim = obs_num_muts_total;
num_trials = 1000;
[ sim_num_muts_per_gene, sim_num_mutsper1kb_per_gene, sim_gene_lengths, sim_gene_genenums, sim_prob_gene_mutated ] = ...
    compute_expected_parallelevo( dir_ref_genome, mut_spec_prob, num_muts_to_sim, num_trials, cdsindex_proteins );





%%

% Histogram of gene lengths
% figure(10); histogram(sim_gene_lengths); xlabel('gene lengths'); ylabel('num genes'); title('all proteins'); set(gca,'FontSize',20)

% if isequal( data_set_name, 'all-muts' )
%     dir_save = dir_results; 
%     plot_file_name = [ 'parallel_' data_set_name '_doublebar-pval-obs-sim' ];
%     plot_compare_pvals_obs_sim( num_genes_mutated, sim_gene_genenums, obs_mutated_genes_genenums, ...
%         obs_num_muts_per_gene, num_proteins_in_genome, num_muts_to_sim, num_trials, sim_num_muts_per_gene, sim_prob_gene_mutated, dir_save, plot_file_name )
% end


%% COMPARE OBSERVED VS EXPECTED NUMBERS OF MULTIPLY MUTATED GENES
fprintf(1,['Computing dN/dS as a function of mutational density...' '\n'])





% Filters
min_num_muts_per_gene = 2; % only do density for multiply mutated genes
list_density_thresholds = 1:1:5; 

% Number of gene sets to investigate
num_genesets_of_interest = numel(list_density_thresholds)+1;

% % Initialize
genesets_of_interest_genenums = cell( numel(list_density_thresholds)+1,1 );

% Find genes of interest
for i=1:numel(list_density_thresholds)
    this_density_threshold = list_density_thresholds(i);
    genes_of_interest_bool = ( ...
        obs_num_muts_per_gene >= min_num_muts_per_gene ...
        & obs_num_mutsper1kb_per_gene >= this_density_threshold ...
        );
    genesets_of_interest_genenums{i} = obs_mutated_genes_genenums( genes_of_interest_bool );
end
genesets_of_interest_genenums{end} = obs_mutated_genes_genenums; % all genes

% Get dN/dS info for genes of interest
[ ~, ~, dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset, ~ ] = ...
    get_dNdS_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_geneinfo );

[dstop_by_geneset, dstoplower_by_geneset, dstop_upper_by_geneset] = ...
    get_dstop_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_geneinfo );


%% Figure: dN/dS for genes of interest



% 
% 
inf_substitute_for_plotting = 100; % since we can't plot infinity
dNdS_for_plot = dNdS_by_geneset;
dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ];
CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
plot_title = 'dN/dS by mutational density';
x_axis_label = 'mutational density threshold (mut/kb)';
density_list_labels = arrayfun(@(x) { num2str(x) }, list_density_thresholds );
density_list_labels{end+1} = 'all';
x_tick_labels = density_list_labels;
dir_save = dir_results ; 
plot_file_name = [ 'parallel-sets_' data_set_name '_dNdS_density' ];
y_axis_factor = 4;
% Call plot function
plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 

% Save data
save( [ dir_results '/' 'data_dNdS_by_density.mat' ], ...
    'min_num_muts_per_gene', 'list_density_thresholds', ...
    'dNdS_by_geneset', 'dNdS_lower_by_geneset', 'dNdS_upper_by_geneset' )

%%




%%






% inf_substitute_for_plotting = 100; % since we can't plot infinity
% dSTOP_for_plot = dstop_by_geneset;
% dSTOP_for_plot( dSTOP_for_plot==Inf ) = inf_substitute_for_plotting;
% CI_for_plot = [ dstoplower_by_geneset; dstop_upper_by_geneset ];
% CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
% plot_title = 'dSTOP/dNONSTOP by mutational density';
% x_axis_label = 'mutational density threshold (mut/kb)';
% density_list_labels = arrayfun(@(x) { num2str(x) }, list_density_thresholds );
% density_list_labels{end+1} = 'all';
% x_tick_labels = density_list_labels;
% dir_save = dir_results ; 
% plot_file_name = [ 'parallel-sets_' data_set_name '_dSTOP_density' ];
% y_axis_factor = 20;
% % Call plot function
% plot_dnds_bar( dSTOP_for_plot, CI_for_plot, ...
%     plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
%     dir_save, plot_file_name ) 
% 





%% IDENTIFY GENES OF INTEREST USING POISSON STATISTICS AND BENJAMINI-HOCHBERG PROCEDURE





% Poisson pvalue for the observed number of mutations on each gene
pval_poisson = zeros( num_genes_mutated,1 );
for n=1:num_genes_mutated
    pval_poisson(n) = poisscdf( obs_num_muts_per_gene(n)-1, ...
        sim_prob_gene_mutated( obs_mutated_genes_genenums(n)==sim_gene_genenums )*num_muts_to_sim, ...
        'upper' ); % poisscdf( x, labmda, 'upper' ) computes pval for >x not >=x
end

bh_thresholds = benjamini_hochberg_simple( pval_poisson, num_proteins_in_genome, fdr_to_test );





%% Make a big table summarizing all genes mutated





annotations_mutgenes = struct;

for g=1:numel(obs_mutated_genes_genenums)
    % Basic info about gene and its mutations
    genecluster_mut_idx = find([annotation_full_everything.gene_num]==obs_mutated_genes_genenums(g));
    annotations_mutgenes(g).lineage = str2double(extractAfter(data_set_name,'_'));
    annotations_mutgenes(g).genenum = obs_mutated_genes_genenums(g);
    annotations_mutgenes(g).num_muts_obs = obs_num_muts_per_gene(g);
    % Poisson pvalue
    annotations_mutgenes(g).pval_poisson = pval_poisson(g);
    % Benjamini-Hochberg value
    annotations_mutgenes(g).BH_threshold = bh_thresholds(g);
    % Observed N/S
    next_types = [annotation_full_everything(genecluster_mut_idx).type];
    next_num_N = sum(next_types=='N');
    next_num_S = sum(next_types=='S');
    next_num_stop= sum(arrayfun(@(x) endsWith(x{:},'*') , [annotation_full_everything(genecluster_mut_idx).muts]));
    [ p_nonsyn_observed, CI_nonsyn ] = binofit( next_num_N, next_num_N + next_num_S, 0.05 );
    [ p_STOP_observed, CI_STOP ] = binofit( next_num_stop, next_num_N + next_num_S, 0.05 );

    NS_observed = p_nonsyn_observed/(1-p_nonsyn_observed);
    STOP_observed=p_STOP_observed/(1-p_STOP_observed);
    NS_obs_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) );
    NS_obs_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) );
    STOP_obs_lower = ( CI_STOP(1)/(1-CI_STOP(1)) );
    STOP_obs_upper = ( CI_STOP(2)/(1-CI_STOP(2)) );
    annotations_mutgenes(g).num_N = next_num_N;
    annotations_mutgenes(g).num_S = next_num_S;
    annotations_mutgenes(g).num_STOP = next_num_stop;
    annotations_mutgenes(g).NS_obs = NS_observed;
    annotations_mutgenes(g).STOP_observed = STOP_observed;

    annotations_mutgenes(g).NS_obs_lower = NS_obs_lower;
    annotations_mutgenes(g).NS_obs_upper = NS_obs_upper;
    annotations_mutgenes(g).STOP_obs_lower = STOP_obs_lower;
    annotations_mutgenes(g).STOP_obs_upper = STOP_obs_upper;
    % Expected N/S
    p_nonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, obs_mutated_genes_genenums(g) );
    p_stop_expected = compute_probstop( dir_ref_genome, mut_spec_prob, obs_mutated_genes_genenums(g) );

    NS_expected = p_nonsyn_expected./(1-p_nonsyn_expected);
    STOP_expected = p_stop_expected./(1-p_stop_expected);

    annotations_mutgenes(g).NS_expected_gene = NS_expected;
    annotations_mutgenes(g).STOP_expected = STOP_expected;

    annotations_mutgenes(g).NS_expected_lineage = NS_expected_allproteins;
    annotations_mutgenes(g).STOP_expected_allproteins = STOP_expected_allproteins;
    % dN/dS
    annotations_mutgenes(g).dNdS_obs = NS_observed/NS_expected;
    annotations_mutgenes(g).dNdS_obs_lower = NS_obs_lower/NS_expected;
    annotations_mutgenes(g).dNdS_obs_upper = NS_obs_upper/NS_expected;

    annotations_mutgenes(g).dSTOP_observed = STOP_observed/STOP_expected;
    annotations_mutgenes(g).dSTOP_obs_lower = STOP_obs_lower/STOP_expected;
    annotations_mutgenes(g).dSTOP_obs_upper = STOP_obs_upper/STOP_expected;
    % Annotation and protein sequence - just grab the first one
    if numel(genecluster_mut_idx)>1
        next_annotation = annotation_full_everything(genecluster_mut_idx(1)).protein;
        next_nt_sequence = annotation_full_everything(genecluster_mut_idx(1)).Sequence;
        next_aa_sequence = annotation_full_everything(genecluster_mut_idx(1)).translation;
        next_protein_cluster = annotation_full_everything(genecluster_mut_idx(1)).locustag_CDHITcluster;
    else
        next_annotation = annotation_full_everything(genecluster_mut_idx).protein;
        next_nt_sequence = annotation_full_everything(genecluster_mut_idx(1)).Sequence;
        next_aa_sequence = annotation_full_everything(genecluster_mut_idx).translation;
        next_protein_cluster = annotation_full_everything(genecluster_mut_idx).locustag_CDHITcluster;
    end
    annotations_mutgenes(g).protein = next_annotation;
    annotations_mutgenes(g).sequence = next_nt_sequence;
    annotations_mutgenes(g).translation = next_aa_sequence;
    annotations_mutgenes(g).gene_length = obs_genes_lengths(g);
    annotations_mutgenes(g).protein_cluster = next_protein_cluster;
end






%% ANALYSIS

%% For supplemental figure: 
% % number of genes detected for different FDR thresholds
% % dN/dS for top genes

% Non-changing params
min_num_muts_per_gene = 2; % minimum number of mutations on a gene

% Number of genes detected at different FDR's 
num_genes_detected = zeros( size( fdr_to_test ) );

for i=1:numel(fdr_to_test)

    % Parameters
    fdr_bh = fdr_to_test(i); % false discovery rate
    % Strings to record params in fig names
    fdr_str = [ '_fdr-' num2str(100*fdr_bh) ];
    nmuts_str = [ '_minmuts-' num2str(min_num_muts_per_gene) ];
    params_str = [ fdr_str nmuts_str ];
    fprintf(1,['Params: ' params_str '\n'])

    % Benjamini-Hochberg procedure
    plot_file_name = [ 'parallel-BH_' data_set_name '_benjhoch' fdr_str nmuts_str ];
    dir_save = dir_results; 
    [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh, dir_save, plot_file_name );

    % Genes of interest
    genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_per_gene );
    num_genes_of_interest = sum( genes_of_interest_bool );

    num_genes_detected(i) = num_genes_of_interest;
        
end

% dN/dS for top N genes

[ ~, pval_order ] = sort( pval_poisson, 'ascend' );
dNdS_top_genes = zeros( numel(top_num_genes_to_test),1 ); % initialize
CIs_top_genes = zeros( numel(top_num_genes_to_test),2 ); % initialize
genenums_to_save = cell( numel(top_num_genes_to_test),1 );
for n=1:numel(top_num_genes_to_test)
    top_num_genes_to_test_num = top_num_genes_to_test(n);
    top_gene_indices = pval_order( 1:top_num_genes_to_test_num );
    % Genes of interest
    genes_of_interest_genenums = obs_mutated_genes_genenums( top_gene_indices );
    genenums_to_save{n} = genes_of_interest_genenums;
    genes_of_interest_lengths = obs_genes_lengths( top_gene_indices );
    num_genes_of_interest = top_num_genes_to_test_num;
    % Get dN/dS for these genes
    [ ~, ~, dNdS_by_gene, dNdS_lower_by_gene, dNdS_upper_by_gene, ~, ~ ] = ...
        get_dNdS_for_genes_of_interest( num_genes_of_interest, genes_of_interest_genenums, annotation_full_everything, ...
        dir_ref_genome, mut_spec_prob, dir_geneinfo );
    dNdS_top_genes(n) = dNdS_by_gene(end);
    CIs_top_genes(n,:) = [ dNdS_lower_by_gene(end); dNdS_upper_by_gene(end) ];
end
% Save data for supp
save( [ dir_results '/' 'data_for_supp.mat' ], ...
    'fdr_to_test', 'min_num_muts_per_gene', 'top_num_genes_to_test', ...
    'num_genes_detected', 'dNdS_top_genes', 'CIs_top_genes', 'genenums_to_save' )


%% Broad scan through FDRs to identify genes of interest
fprintf(1,['Identifying genes of interest...' '\n'])

fdr_bh_list = fdr_to_test;
min_num_muts_list = [ 2 ];

for i=1:numel(fdr_bh_list)
    for j=1:numel(min_num_muts_list)
        
        % Parameters
        fdr_bh = fdr_bh_list(i); % false discovery rate
        min_num_muts_per_gene = min_num_muts_list(j); % minimum number of mutations on a gene
        % Strings to record params in fig names
        fdr_str = [ '_fdr-' num2str(100*fdr_bh) ];
        nmuts_str = [ '_minmuts-' num2str(min_num_muts_per_gene) ];
        params_str = [ fdr_str nmuts_str ];
        fprintf(1,['Params: ' params_str '\n'])

        % Benjamini-Hochberg procedure
        plot_file_name = [ 'parallel-BH-GOI_' data_set_name '_benjhoch' fdr_str ];
        dir_save = dir_results; 
        [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh, dir_save, plot_file_name );

        % Genes of interest
        genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_per_gene );
        num_genes_of_interest = sum( genes_of_interest_bool );

        genes_of_interest_genenums = obs_mutated_genes_genenums( genes_of_interest_bool );
        genes_of_interest_lengths = obs_genes_lengths( genes_of_interest_bool );

        % Plots describing selection of genes of interest
        % Extra plot: scatter plot of poisson pvalue vs gene length with genes of interest highlighted
        plot_title = 'all within-lineage mutated genes';
        dir_save = dir_results; 
        plot_file_name = [ 'parallel-BH-GOI_' data_set_name '_scatter' params_str ];
        compare_filter_to_poisson( obs_genes_lengths, pval_poisson, obs_mutated_genes_genenums, genes_of_interest_genenums, ...
            plot_title, dir_save, plot_file_name );
        
        if num_genes_of_interest>0

            % Extra plot: histogram of poisson pvalues with genes of interest highlighted
            plot_title = 'all within-lineage mutated genes';
            x_axis_label = 'poisson pvalue';
            y_axis_label = 'num of genes';
            line_label = [ 'pval_cutoff=' num2str(pval_max) ' (BH)' ];
            plot_file_name = [ 'parallel-BH-GOI_' data_set_name '_hist' params_str ];
            plot_parallelevo_hist_exp_line_obs_poisson( ...
                log10(pval_max), log10(pval_poisson), genes_of_interest_bool, ...
                plot_title, x_axis_label, y_axis_label, line_label, dir_results, plot_file_name );

            if num_genes_of_interest<=50

                % Get dN/dS info for genes of interest
                [ Nobs_by_gene, Sobs_by_gene, dNdS_by_gene, dNdS_lower_by_gene, dNdS_upper_by_gene, name_by_gene, pval_dNdS ] = ...
                    get_dNdS_for_genes_of_interest( num_genes_of_interest, genes_of_interest_genenums, annotation_full_everything, ...
                    dir_ref_genome, mut_spec_prob, dir_geneinfo );
                length_by_gene = [ genes_of_interest_lengths, mean(genes_of_interest_lengths) ];

                % Figure: dN/dS for genes of interest
                inf_substitute_for_plotting = 100; % since we can't plot infinity
                dNdS_for_plot = dNdS_by_gene;
                dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
                CI_for_plot = [ dNdS_lower_by_gene; dNdS_upper_by_gene ];
                CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
                plot_title = [ 'genes of interest: ' data_set_name ' (p=' num2str(pval_dNdS) ')' ];
                x_axis_label = 'gene annotation';
                x_tick_labels = name_by_gene;
                dir_save = dir_results ; 
                plot_file_name = [ 'parallel-BH-GOI_' data_set_name '_dNdS' params_str  ];
                text_labels_1 = arrayfun(@(i) [ num2str(Nobs_by_gene(i)) '/' num2str(Sobs_by_gene(i)) ], 1:1:num_genes_of_interest+1, 'UniformOutput', false);
                text_labels_2 = length_by_gene;
                y_axis_factor = 8;
                % Call plot function
                plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
                    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
                    dir_save, plot_file_name, ...
                    text_labels_1, text_labels_2 ) 

            end
            
        end

    end
end


%% Broad scan through FDRs to compute dN/dS
fprintf(1,['Evaluating dN/dS...' '\n'])

fdr_bh_list = fdr_to_test;
min_num_muts_list = [1 2 5 10 ];

% Save dN/dS for all mutations meeting min_num_muts_list criteria
dNdS_return = zeros( size( min_num_muts_list ) );
CIs_return = zeros( numel(min_num_muts_list), 2);

for i=1:numel( min_num_muts_list )
    
    % Identify genes of interest
    num_genesets_of_interest = numel( fdr_bh_list );
    genesets_of_interest_genenums = cell( num_genesets_of_interest,1 );
    genesets_of_interest_size = zeros( num_genesets_of_interest,1 );
    % Loop through BH FDR cutoffs
    for j=1:numel( fdr_bh_list )
        [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh_list(j) );
        genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_list(i) );
        genes_of_interest_genenums = obs_mutated_genes_genenums( genes_of_interest_bool );
        genesets_of_interest_genenums{j} = genes_of_interest_genenums;
        genesets_of_interest_size(j) = sum( genes_of_interest_bool );
    end

    % Get dN/dS info for genes of interest
    [ Nobs_by_geneset, Sobs_by_geneset, dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset, dnds_by_geneset ] = ...
        get_dNdS_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
        dir_ref_genome, mut_spec_prob, dir_geneinfo );

    % Figure: dN/dS for genes of interest
    inf_substitute_for_plotting = 100; % since we can't plot infinity
    dNdS_for_plot = dNdS_by_geneset;
    dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
    CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ];
    CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
    plot_title = 'genes of interest: dN/dS';
    x_axis_label = 'BH FDR';
    fdr_bh_list_labels = arrayfun(@(x) { num2str(x) }, fdr_bh_list );
    fdr_bh_list_labels{end} = 'all';
    x_tick_labels = fdr_bh_list_labels;
    dir_save = dir_results ; 
    plot_file_name = [ 'parallel-sets_' data_set_name '_dNdS_nmuts-' num2str(min_num_muts_list(i)) ];
    text_labels_1 = dnds_by_geneset;
    text_labels_2 = genesets_of_interest_size;
    y_axis_factor = 4;
    % Call plot function
    plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name, ...
        text_labels_1, text_labels_2 ) 
    
    % Save
    dNdS_return(i) = dNdS_for_plot(end);
    CIs_return(i,:) = CI_for_plot(:,end);
    
end



end
