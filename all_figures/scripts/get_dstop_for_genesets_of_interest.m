function [dSTOP_by_geneset, dSTOP_lower_by_geneset, dSTOP_upper_by_geneset ] = ...
    get_dstop_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_save ) 


%% Initialize

dSTOP_by_geneset  = zeros( 1,num_genesets_of_interest ); % initialize
dSTOP_lower_by_geneset = zeros( 1,num_genesets_of_interest ); % initialize
dSTOP_upper_by_geneset = zeros( 1,num_genesets_of_interest ); % initialize


%% dN/dS for potential genes of interest

% Loop through genes of interest
for n=1:numel(genesets_of_interest_genenums)
    
    % Get gene number
    my_genenums = genesets_of_interest_genenums{n};
    
    % Compute observed N/S
    [ p_STOP, CI_STOP, num_muts_NONSTOP, num_muts_STOP ] = ...
        compute_observed_dSTOP( annotation_full_everything, my_genenums );
    dSTOP_observed = (p_STOP/(1-p_STOP)); % N/S observed
    fprintf(1,[ 'NONSTOP: ' num2str(num_muts_NONSTOP) '\n' ])
    fprintf(1,[ 'STOP: ' num2str(num_muts_STOP) '\n' ])
    
    % Compute expected N/S
    probSTOP_expected = compute_probstop( dir_ref_genome, mut_spec_prob, my_genenums );
    dSTOP_expected_local = probSTOP_expected/(1-probSTOP_expected); % N/S expected from neutral model
    
    % Compute dN/dS
    dSTOP = dSTOP_observed/dSTOP_expected_local; % relative to neutral model
    CI_lower = ( CI_STOP(1)/(1-CI_STOP(1)) ) / dSTOP_expected_local;
    CI_upper = ( CI_STOP(2)/(1-CI_STOP(2)) ) / dSTOP_expected_local;
    
    % Record
    dSTOP_by_geneset(n)  = dSTOP;
    dSTOP_lower_by_geneset(n) = CI_lower;
    dSTOP_upper_by_geneset(n) = CI_upper;
    
    fprintf(1,'\n')
    
end


end