function all_proteins_info = get_all_protein_info(uclades,mut_spec_prob,dir_genomes)
all_proteins_info = struct;
all_proteins_tally = 0;

for clade_idx = 1:numel(uclades)
    clade_number = uclades(clade_idx);
    fprintf([ 'Pulling CDS for Clade ' num2str(clade_number) '...\n' ])

    % Where to find reference genomes
    dir_ref_genome = [ dir_genomes num2str(clade_number) ];


    %% Get gene numbers for proteins only (exclude tRNA and rRNA)
    fprintf(1,['Getting genome information...' '\n'])

    [ proteins_cdsindex, proteins_locustag ] = get_coding_genenums( dir_ref_genome); % avoid genes in regions masked in SNP calling


    %% Expected probability that each gene is mutated
    fprintf(1,['Simulating mutated genes...' '\n'])
[ ~, ~, sim_gene_lengths, ~, sim_prob_gene_mutated ] = compute_expected_parallelevo( dir_ref_genome, mut_spec_prob, proteins_cdsindex );

    
    %% Put protein info into big table

    for p=1:numel(proteins_cdsindex)
        all_proteins_tally = all_proteins_tally+1;
        all_proteins_info(all_proteins_tally).clade = clade_number;
        all_proteins_info(all_proteins_tally).cds_idx = proteins_cdsindex(p);
        all_proteins_info(all_proteins_tally).locustag = proteins_locustag{p};
        all_proteins_info(all_proteins_tally).gene_length = sim_gene_lengths(p);
        all_proteins_info(all_proteins_tally).prob_gene_mutated = sim_prob_gene_mutated(p);
    end

end

end