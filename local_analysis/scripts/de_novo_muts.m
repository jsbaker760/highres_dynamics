function [LineageData,cacnes_annotations_all,sepi_annotations_all]=de_novo_muts(SamplesCSV,dir_lineage_coassemblies,dir_phylogenies)

%CSV used on cluster to align samples to clade co-assemblies
clade_names = string(unique(SamplesCSV.ReferenceGenome));
% create parameters structure
params = struct;
% % For filtering individual calls
params.min_qual_for_call = 30;
params.min_maf_for_call = .75;
params.min_cov_for_call = 4;
% For filtering positions
params.max_fraction_ambigious_samples = .4;
% to filter maNT for anlysis 
params.min_maf_for_analysis = .7;
params.min_coverage_for_analysis = 5; % for each call (fwd+rev)
% Promoter mutations: how far upstream of the nearest gene to annotate something a promoter mutation
params.promotersize = 150;
params.NTs = 'ATCG';    
% location of asssemblies for each clade with filtered contigs
params.REFGENOMEDIRassemblies = [pwd dir_lineage_coassemblies];
% where to put subject-seperated phylogenies
% directory in which to place S. epi assembly trees


% % identify de novo mutations  in Sepi, build trees, get data for other thing
SepiLineages=identify_mutations_lineage_assemblies(dir_phylogenies,clade_names(contains(clade_names,"sepi")), params);

% identify de novo mutations  in Cacnes, build trees, get data for other thing
CacnesLineages = identify_mutations_lineage_assemblies(dir_phylogenies,clade_names(contains(clade_names,"cacnes")), params);

% get lineage data and put it all into a single structure
SepiLineages = struct2table(SepiLineages);SepiLineages.SpeciesName = repmat("sepi",numel(SepiLineages.cladenumber),1);
CacnesLineages = struct2table(CacnesLineages);CacnesLineages.SpeciesName = repmat("cacnes",numel(CacnesLineages.cladenumber),1);
LineageData = [CacnesLineages ;    SepiLineages ];

% get annotations for both species
sepi_annotations_all = horzcat(LineageData.annotations{LineageData.SpeciesName=="sepi"});
cacnes_annotations_all=horzcat(LineageData.annotations{LineageData.SpeciesName=="cacnes"});