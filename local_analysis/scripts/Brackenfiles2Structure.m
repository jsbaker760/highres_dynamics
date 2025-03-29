%% dir is the directory where bracken files are with the format:
% name	 taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
function BrackenDataAll =  Brackenfiles2Structure(SampleNames,dir)
% dir = 'data/bracken/metagenomics_tubes/';
% for reasons i don't remember, the bracken files in the dir end in .dat
species_suffix = '.S.bracken';
genus_suffix = '.G.bracken';
%% Brackenfiles: a string list of the pertinant filenames in the directory above
% BrackenSampleNames: the name of the files in Brackenfiles without 'BrackenS/2-kraken2/'


species_Brackenfiles = string(strsplit(ls([dir '*' species_suffix])))';
species_Brackenfiles = species_Brackenfiles(~arrayfun(@(x) isempty(x{:}),species_Brackenfiles));
genus_Brackenfiles = string(strsplit(ls([dir '*' genus_suffix])))';
genus_Brackenfiles = genus_Brackenfiles(~arrayfun(@(x) isempty(x{:}),genus_Brackenfiles));

species_BrackenSampleNames = strrep(species_Brackenfiles,dir,'');
species_BrackenSampleNames = strrep(species_BrackenSampleNames,species_suffix,'');
genus_BrackenSampleNames = strrep(genus_Brackenfiles,dir,'');
genus_BrackenSampleNames = strrep(genus_BrackenSampleNames,genus_suffix,'');


%% samples.mat contains variable 'samples', a (#Samplesx1) string array of sample names, usually from samples.csv
% the final data structure BrackenDataAll has the order below
% such that the bracken directory can contain files which are not in BrackenDataAll
% if its name is not in 'samples'

%% BrackenDataAll contains all the bracken data
% each sample is a table which can be accessed with dot indexing i.e.
% Sample1Table = BrackenDataAll(1)
% Sample1Species = BrackenDataAll(1).name
% The reason the structures are numbered instead of
% using samplenames is that in MATLAB fields cannot start with a number
% whereas my samples do

BrackenDataAll = struct;
BrackenDataAll.Species = struct;
BrackenDataAll.Genus = struct;

%% Iterate through samples
for s = 1:numel(SampleNames)
    % find bracken file corrosponding to sample
    has_species_file = SampleNames(s)==species_BrackenSampleNames;
    if sum(has_species_file)==1
        file = species_Brackenfiles{has_species_file};
    	% open file
        fid = fopen(file);
        % parse file into cell with textscan
        t = textscan(fid,'%s%f%s%f%f%f%f','Delimiter','\t','Headerlines',1);
        % close file
        fclose(fid);

        BrackenDataAll.Species(s).name                    =string(t{1});
        BrackenDataAll.Species(s).taxonomy_id             =single(t{2});
        BrackenDataAll.Species(s).taxonomy_lvl            =string(t{3});
        BrackenDataAll.Species(s).kraken_assigned_reads   =single(t{4});
        BrackenDataAll.Species(s).added_reads             =single(t{5});
        BrackenDataAll.Species(s).new_est_reads           =single(t{6});
        BrackenDataAll.Species(s).fraction_total_reads    =single(t{7});
    end
    has_genus_file = SampleNames(s)==genus_BrackenSampleNames;
    if sum(has_genus_file)==1
        file = genus_Brackenfiles{has_genus_file};
    	% open file
        fid = fopen(file);
        % parse file into cell with textscan
        t = textscan(fid,'%s%f%s%f%f%f%f','Delimiter','\t','Headerlines',1);
        % close file
        fclose(fid);

        BrackenDataAll.Genus(s).name                    =string(t{1});
        BrackenDataAll.Genus(s).taxonomy_id             =single(t{2});
        BrackenDataAll.Genus(s).taxonomy_lvl            =string(t{3});
        BrackenDataAll.Genus(s).kraken_assigned_reads   =single(t{4});
        BrackenDataAll.Genus(s).added_reads             =single(t{5});
        BrackenDataAll.Genus(s).new_est_reads           =single(t{6});
        BrackenDataAll.Genus(s).fraction_total_reads    =single(t{7});
    end

end

AllTaxaSpecies = unique(vertcat(BrackenDataAll.Species.name));
AllTaxaGenera  = unique(vertcat(BrackenDataAll.Genus.name));

%%


BrackenDataAll.SpeciesArray = struct;
[BrackenDataAll.SpeciesArray.kraken_assigned_reads,    ...
    BrackenDataAll.SpeciesArray.added_reads,              ...
    BrackenDataAll.SpeciesArray.new_est_reads,            ...
    BrackenDataAll.SpeciesArray.fraction_total_reads        ] = deal(zeros(numel(AllTaxaSpecies),numel(SampleNames)));
%
BrackenDataAll.GenusArray   = struct;
[BrackenDataAll.GenusArray.kraken_assigned_reads,      ...
    BrackenDataAll.GenusArray.added_reads,                ...
    BrackenDataAll.GenusArray.new_est_reads,              ...
    BrackenDataAll.GenusArray.fraction_total_reads          ] = deal(zeros(numel(AllTaxaGenera),numel(SampleNames)));


%%

for s =1:numel(SampleNames)

    idxS=arrayfun(@(x) find(x==AllTaxaSpecies), BrackenDataAll.Species(s).name);
    BrackenDataAll.SpeciesArray.kraken_assigned_reads(idxS,s) = BrackenDataAll.Species(s).kraken_assigned_reads;
    BrackenDataAll.SpeciesArray.added_reads(idxS,s)           = BrackenDataAll.Species(s).added_reads;
    BrackenDataAll.SpeciesArray.new_est_reads(idxS,s)         = BrackenDataAll.Species(s).new_est_reads;
    BrackenDataAll.SpeciesArray.fraction_total_reads(idxS,s)  = BrackenDataAll.Species(s).fraction_total_reads;

    idxG=arrayfun(@(x) find(x==AllTaxaGenera),  BrackenDataAll.Genus(s).name);
    BrackenDataAll.GenusArray.kraken_assigned_reads(idxG,s) = BrackenDataAll.Genus(s).kraken_assigned_reads;
    BrackenDataAll.GenusArray.added_reads(idxG,s)           = BrackenDataAll.Genus(s).added_reads;
    BrackenDataAll.GenusArray.new_est_reads(idxG,s)         = BrackenDataAll.Genus(s).new_est_reads;
    BrackenDataAll.GenusArray.fraction_total_reads(idxG,s)  = BrackenDataAll.Genus(s).fraction_total_reads;
end


%%
BrackenDataAll.SpeciesArray.TaxaNames = AllTaxaSpecies;
BrackenDataAll.GenusArray.TaxaNames   = AllTaxaGenera;

%%

end