function [ proteins_cdsindex, proteins_locustag ] = get_coding_genenums( dir_ref_genome, positions_to_mask )

%% Get CDS
if ~exist([dir_ref_genome '/cds_sorted.mat'], 'file')
    read_gb_assembly(dir_ref_genome)
end % makes cds_sorted.mat if it doesn't already exist
load([dir_ref_genome '/cds_sorted.mat'],'CDS') 
all_genes=[CDS{:}]; % copies CDS






%% Get CDS indices for proteins only
proteins_cdsindex = find([all_genes.tagnumber]>0.5); % excludes tRNA and rRNA
proteins_locustag = { all_genes(proteins_cdsindex).locustag }; 





%% Remove genes containing masked positions
if nargin>1
    [ChrStarts, GLength, ~, ~]= genomestats(dir_ref_genome);
    pos_to_mask = p2chrpos( positions_to_mask, ChrStarts );
    locustags_mask = {};
    for p=1:numel(positions_to_mask)
        CDS_Chr = CDS{pos_to_mask(p,1)};
        if ~isempty(CDS_Chr)
            gene_loc1 = [CDS_Chr.loc1];
            gene_loc2 = [CDS_Chr.loc2];
            cds_chr_index = find( pos_to_mask(p,2)>=gene_loc1 & pos_to_mask(p,2)<=gene_loc2 );
            if ~isempty( cds_chr_index )
                locustags_mask{end+1} = CDS_Chr(cds_chr_index).locustag;
            end
        end
    end
    locustags_mask = unique(locustags_mask);
    proteins_cdsindex = find( ~ismember( proteins_locustag,locustags_mask ) );
    fprintf(1,['Masking ' num2str(numel(proteins_locustag)-numel(proteins_cdsindex)) ' of ' num2str(numel(proteins_locustag)) ' proteins in genome.\n'  ])
end





%% Grab locus tags for unmasked proteins





proteins_locustag = { all_genes(proteins_cdsindex).locustag };





end