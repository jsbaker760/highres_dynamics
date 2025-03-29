function saveSupplementaryTablesExcel(inputData)
% saveSupplementaryTablesExcel saves output data into 6 Excel files (Table S1-S6)
% with multiple sheets instead of separate files.

% If inputData is a file path, load the MAT file; otherwise assume a structure.
if ischar(inputData) || isstring(inputData)
    S = load(inputData);
else
    S = inputData;
end

% Unpack S into workspace variables.
fields = fieldnames(S);
for k = 1:length(fields)
    eval([fields{k} ' = S.(fields{k});']);
end

% Create output folder if it doesn't exist.
if ~exist('SupplementaryTables', 'dir')
    mkdir('SupplementaryTables');
end

%% (Preliminary data processing)

% Table S2: sheet cdhit_homologues_isolates
SampleNames = CDHITData.SampleNames';
NumClusters = size(CDHITData.CDHITClusterxSamples, 1);
ClusterNumbers = 0:NumClusters-1; % cd-hit is zero indexed
Tsamples = arrayfun(@(x) {strjoin(string(ClusterNumbers(CDHITData.CDHITClusterxSamples(:,x))), " ")}, 1:numel(SampleNames));
CDHitClusterTable = table;
CDHitClusterTable.SampleNames = SampleNames';
CDHitClusterTable.CDHIT_clusters_found_in_sample = vertcat(Tsamples{:});

% Table S3: sheet abundance_subject_timepoint_loc
MtubesTable = parse_to_csv_version(Mtubes);
% Table S3: sheet abundance_subject_timepoint
MSubjectTimeTable = parse_to_csv_version(MsubjectTime);
% Table S3: sheet abundance_subject
MSubjectTable = parse_to_csv_version(Msubject);
% parse illumina barcodes
MetagenomicsSampleNames.PROVIDER_ls=string(MetagenomicsSampleNames.PROVIDER_ls);
BCs_tubes=arrayfun(@(x) {string(MetagenomicsSampleNames.PROVIDER_ls(ismember(MetagenomicsSampleNames.SAMPLE_ls,x{:})))}, Mtubes.SampleNames);
BCs_plates=arrayfun(@(x) {string(MetagenomicsSampleNames.Plate(ismember(MetagenomicsSampleNames.SAMPLE_ls,x{:})))}, Mtubes.SampleNames);
BCs_wells=arrayfun(@(x) {string(MetagenomicsSampleNames.Well(ismember(MetagenomicsSampleNames.Index,x{:})))}, Mtubes.index604);
T_barcodes = table;
Ndiff = cellfun(@numel,BCs_wells);
for i = 1:numel(BCs_tubes)
    N = Ndiff(i);
    T_barcodes.SequencingPlate(i) = strjoin(BCs_plates{i},'\t');
    T_barcodes.SequencingWell(i) = strjoin(BCs_wells{i},'\t');
    T_barcodes.illuminaBarcodes(i) = strjoin(BCs_tubes{i},'\t');
end
% Build the TaxonNames table using the first element of each variable
taxonNames = table;
taxonNames.CacnesPhylotypeNames   = {strjoin(string(MsubjectTime.CacnesPhylotypeNames(1,:)), '\t')};
taxonNames.SepiPhylotypeNames      = {strjoin(string(1:4), '\t')};
taxonNames.SepiSubPhylotypeNames   = {strjoin(string(MsubjectTime.SepiSubphylotypeNames(1,:)), '\t')};
taxonNames.SepiLineageNumbers      = {strjoin(string(MsubjectTime.SepiLineageNumbers(1,:)), '\t')};
taxonNames.CacnesLineageNumbers    = {strjoin(string(MsubjectTime.CacnesLineageNumbers(1,:)), '\t')};
taxonNames.CacnesSubphylotypesNames= {strjoin(["A_1", "A_2", "A_3", "B_1", "C_1", "D_1", "E_1", "F_1", "F_2", "F_3", "H_1", "H_2", "H_3", "K_1_1", "K_1_2", "K_2", "K_3", "L_1"], '\t')};
taxonNames.BrackenSpeciesNames     = {strjoin(string(MsubjectTime.NamesSpecies(1,:)), '\t')};
taxonNames.BrackenGenusNames       = {strjoin(string(MsubjectTime.NamesGenus(1,:)), '\t')};
taxonNames.BrackenSpeciesNamesCutotypeClustering  = {strjoin(string(MsubjectTime.BrackenNamesClustering(1,:)), '\t')};
% parse variable order
MtubesTable=[MtubesTable T_barcodes];
% split up files to not print too many names
MtubesTable = MtubesTable(:,~ismember(string(MtubesTable.Properties.VariableNames),["CacnesPhylotypeNames","SepiPhylotypeNames","SepiSubPhylotypeNames","SepiLineageNumbers","CacnesLineageNumbers","SpeciesNamesBracken","GenusNamesBracken","BrackenNamesClustering"]));
MSubjectTimeTable = MSubjectTimeTable(:,~ismember(string(MSubjectTimeTable.Properties.VariableNames),["CacnesPhylotypeNames","SepiPhylotypeNames","SepiSubPhylotypeNames","SepiLineageNumbers","CacnesLineageNumbers","SpeciesNamesBracken","GenusNamesBracken","BrackenNamesClustering"]));
MSubjectTable = MSubjectTable(:,~ismember(string(MSubjectTable.Properties.VariableNames),["CacnesPhylotypeNames","SepiPhylotypeNames","SepiSubPhylotypeNames","SepiLineageNumbers","CacnesLineageNumbers","SpeciesNamesBracken","GenusNamesBracken","BrackenNamesClustering"]));

% Table S4 coassembly statistics
T_coassembly_stats=table;
% get clades which were used in analysis and not excluded so other tables
% can print them in the same order
for i = 1:176
    CladeName = Clades(i);
    F = ['data/clade_coassemblies/' char(CladeName) '/genome.fasta'];
    sequences=fastaread(F);
    Headers =struct2table(sequences);
    Headers=string(Headers.Header);
    Headers=arrayfun(@(x) {strsplit(x,"cov_")} , Headers);
    Headers=arrayfun(@(x) str2double(x{:}(2)) , Headers);
    seqLengths = cellfun(@numel, {sequences.Sequence});
    sortedLengths = sort(seqLengths, 'descend');
    totalLength = sum(sortedLengths);
    % Compute the N50
    cumulativeSum = cumsum(sortedLengths);
    n50 = sortedLengths(find(cumulativeSum >= totalLength / 2, 1));
    T_coassembly_stats.CladeName(i)=CladeName;
    T_coassembly_stats.CoverageCutoff(i) = CutoffHeights(i);
    T_coassembly_stats.N50(i)=n50;
    T_coassembly_stats.Length(i)=totalLength;
    T_coassembly_stats.MinCoverageAnyContig(i)=min(Headers);
end
% sort them numerically
parsednames = arrayfun(@(x) {strsplit(x,"_")}, T_coassembly_stats.CladeName);
parsednames = vertcat(parsednames{:});
[~,i]=sort(str2double(parsednames(:,3)));
parsednames=parsednames(i,:);
T_coassembly_stats=T_coassembly_stats(i,:);
[~,i]=sort(parsednames(:,1));
T_coassembly_stats=T_coassembly_stats(i,:);

% sheet clade_all_nut_positions
LineageNTsSupplement = table;
LineageNTsSupplement.AncAlleles= arrayfun(@(x) {LineageData.anc_nti{x}(LineageData.goodpos{x})}, 1:size(LineageData,1))';
LineageNTsSupplement.Samplenames= LineageData.samplenames;
LineageNTsSupplement.in_outgroup= LineageData.outgroup;
LineageNTsSupplement.p= LineageData.p;
LineageNTsSupplement.alleles= LineageData.calls_analysis;
LineageNTsSupplement.goodpos = LineageData.goodpos;
R = size(LineageNTsSupplement,1);
NTs=["N" "A" "T" "C" "G"];
names_clades="";
tbs_nucleotides = cell(R,1);
for i = 1:R
    row = LineageNTsSupplement(i,:);
    name = strjoin([LineageData.SpeciesName(i) "clade" string(LineageData.cladenumber(i))],'_');
    names_clades(i)=name;
    newtable = table;
    newtable.name = string(row.Samplenames{:})';
    newtable.is_outgroup_sample = row.in_outgroup{:}';
    nucleotide = row.alleles{:};
    nucleotide=NTs(nucleotide+1)';
    pos = row.p{:}(row.goodpos{:}'==1)';
    for j = 1:numel(pos)
        newtable.(strjoin(["nt_position_" string(pos(j))],''))=nucleotide(:,j);
    end
    newtable.name(end+1)="ANSCESTRAL_ALLELE";
    anc_alleles = row.AncAlleles{:}';
    anc_alleles=NTs(anc_alleles+1);
    try
        newtable(end,3:end)=cellstr(anc_alleles);
    catch
    end
    tbs_nucleotides{i}=newtable;
end
% sort to keep in same order
parsednames = arrayfun(@(x) {strsplit(x,"_")}, names_clades);
parsednames = vertcat(parsednames{:});
[~,i]=sort(str2double(parsednames(:,3)));
parsednames=parsednames(i,:);
names_clades=names_clades(i);
tbs_nucleotides=tbs_nucleotides(i);
[~,i]=sort(parsednames(:,1));
names_clades=names_clades(i);
tbs_nucleotides=tbs_nucleotides(i);

% sheet all_gain_loss_regions
gain_loss_all =cell(R,1);
for i =1:numel(names_clades)
    GainLossRegions = ['data/MGEs/Trees/' char(names_clades(i)) '/' char(names_clades(i)) '.csv'];
    tbl = ['data/MGEs/Regions/' char(names_clades(i)) '/' char(names_clades(i)) '_regions.mat'];
    if isfile(tbl)&isfile(GainLossRegions)
        T = load(tbl);
        T=struct2table(T.all_regions_table);
        T=parse_gain_loss(T);
        gain_loss_all{i}=(T);
    end
end
gain_loss_all(cellfun(@isempty, gain_loss_all))=[];
gain_loss_all=vertcat(gain_loss_all{:});

%% TABLE S1
% Sheet: subject timepoint metadata (SampleInfo)
fileS1 = fullfile('SupplementaryTables', 'Table_S1.xlsx');
writetable(SampleInfo, fileS1, 'Sheet', 'subject_timepoint_metadata');

%% TABLE S2
% Sheets:
%   ncbi_reference_genomes: ncbi_references
%   assembled_isolate_genomes: AssemblyStats
%   cdhit_homologues_isolates: CDHit data
fileS2 = fullfile('SupplementaryTables', 'Table_S2.xlsx');
writetable(ncbi_references, fileS2, 'Sheet', 'ncbi_reference_genomes');
writetable(AssemblyStats([1 3:end],:), fileS2, 'Sheet', 'assembled_isolate_genomes'); %#ok<*USENS>
writetable(CDHitClusterTable, fileS2, 'Sheet', 'cdhit_homologues_isolates');

%% TABLE S3
% Sheets:
%   abundance_subject_timepoint_loc: MtubesTable
%   abundance_subject_timepoint: MSubjectTimeTable
%   abundance_subject: MSubjectTable
%   TaxonNames: tab‐delimited strings for a set of taxon name variables (from MsubjectTime)
fileS3 = fullfile('SupplementaryTables', 'Table_S3.xlsx');
writetable(MtubesTable, fileS3, 'Sheet', 'abundance_subject_timepoint_loc');
writetable(MSubjectTimeTable, fileS3, 'Sheet', 'abundance_subject_timepoint');
writetable(MSubjectTable, fileS3, 'Sheet', 'abundance_subject');
writetable(taxonNames, fileS3, 'Sheet', 'TaxonNames');

%% TABLE S4
% Sheets:
%   clade_coassembly_stats: loaded below
%   lineage_membership: LineageData
%   subject_lineage_dmrcas: LineageDMRCAsForSupplement
%   all_possible_transmissions: AllPossibleTransmissionPatterns
%   all_gain_loss_regions: gain_loss_all
%   clade_all_nut_positions: nucleotide positions per clade
fileS4 = fullfile('SupplementaryTables', 'Table_S4.xlsx');

% sheet clade_coassembly_stats
writetable(T_coassembly_stats, fileS4, 'Sheet', 'clade_coassembly_stats');

% sheet lineage_membership
SIDs=arrayfun(@(x) strjoin(LineageDMRCAs.SID(LineageDMRCAs.Species==LineageData.SpeciesName(x)&LineageDMRCAs.LineageNumber==LineageData.cladenumber(x)),' '), 1:size(LineageData,1)); %#ok<*NODEF>
Nisolates=arrayfun(@(x) strjoin(string(LineageDMRCAs.Nisolates(LineageDMRCAs.Species==LineageData.SpeciesName(x)&LineageDMRCAs.LineageNumber==LineageData.cladenumber(x))),' '), 1:size(LineageData,1));
LineageData.Subjects_with_isolates=SIDs';
LineageData.isolates_per_subjects=Nisolates';
writetable(LineageData(:,15:19), fileS4, 'Sheet', 'lineage_membership');

% sheet subject_lineage_mrcas
LineageSubjectDMRCAsForSupplement = LineageData2supplement(LineageDMRCAs);
writetable(LineageSubjectDMRCAsForSupplement, fileS4, 'Sheet', 'subject_lineage_dmrcas');

% sheet all_possible_transmissions
AllPossibleTransmissionPatternsCacnes.species = repmat("cacnes",size(AllPossibleTransmissionPatternsCacnes,1),1);
AllPossibleTransmissionPatternsSepi.species = repmat("sepi",size(AllPossibleTransmissionPatternsSepi,1),1);
AllPossibleTransmissionPatterns = [AllPossibleTransmissionPatternsCacnes; AllPossibleTransmissionPatternsSepi];
AllPossibleTransmissionPatterns.Recipients=arrayfun(@(x) strjoin(x{:}," "), AllPossibleTransmissionPatterns.Recipients);
AllPossibleTransmissionPatterns.Ntransmissions=arrayfun(@(x) strjoin(string(x{:})," "), AllPossibleTransmissionPatterns.Ntransmissions);
AllPossibleTransmissionPatterns.NisolatesRecipient=arrayfun(@(x) strjoin(string(x{:})," "), AllPossibleTransmissionPatterns.NisolatesRecipient);
AllPossibleTransmissionPatterns.NisolatesSource=arrayfun(@(x) strjoin(string(x{:})," "), AllPossibleTransmissionPatterns.NisolatesSource);
writetable(AllPossibleTransmissionPatterns, fileS4, 'Sheet', 'all_possible_transmissions');

% sheet all_gain_loss_regions
writetable(gain_loss_all, fileS4, 'Sheet', 'all_gain_loss_regions');

% write nucleotide data
% Determine the maximum number of columns among all clade tables.
maxWidth = max(cellfun(@(T) width(T), tbs_nucleotides));
% Initialize a cell array to hold the combined data.
combinedData = {};
for i = 1:numel(names_clades)
    % Create a header row for the clade name. 
    % Place the clade name in column A and pad with empty strings.
    cladeHeader = repmat({''}, 1, maxWidth);
    cladeHeader{1} = names_clades{i};
    combinedData = [combinedData; cladeHeader];
    % Get the table header (variable names) and pad to maxWidth.
    T = tbs_nucleotides{i};
    tblHeader = T.Properties.VariableNames;
    if numel(tblHeader) < maxWidth
        tblHeader = [tblHeader, repmat({''}, 1, maxWidth - numel(tblHeader))];
    end
    combinedData = [combinedData; tblHeader];
    % Convert the table to cell array and pad each row.
    tblData = table2cell(T);
    for j = 1:size(tblData, 1)
        rowData = tblData(j, :);
        if numel(rowData) < maxWidth
            rowData = [rowData, repmat({''}, 1, maxWidth - numel(rowData))];
        end
        combinedData = [combinedData; rowData];
    end
    % Optionally, add an empty row as a separator between clades.
    combinedData = [combinedData; repmat({''}, 1, maxWidth)];
end
% Write the combined data to the Excel sheet.
writecell(combinedData, fileS4, 'Sheet', 'clade_all_nut_positions');

%% TABLE S5
% Sheets:
%   Identity_of_each_clustered_iso: table IT
%   samples_with_outgroups: T_all(:,[2 3 5])
fileS5 = fullfile('SupplementaryTables', 'Table_S5.xlsx');
writetable(IT, fileS5, 'Sheet', 'Identity_of_each_clustered_iso');
writetable(T_all(:, [2:5]), fileS5, 'Sheet', 'samples_with_outgroups');

%% TABLE S6
% Sheets:
%   cacnes_annotations_all: cacnes_annotations_all
%   semi_annotations_all: sepi_annotations_all
%   HomologMembershipCacnes: HomologousClustersCacnesAssemblies
%   HomologMembershipSepi: HomologousClustersSepiAssemblies
%   CacnesAllProteinInfo: CacnesAllProteinInfo
%   SepiAllProteinInfo: SepiAllProteinInfo
%   CacnesAnnotationsHomologs: CacnesAnnotationsHomologues
%   SepiAnnotationsHomologs: SepiAnnotationsHomologues
fileS6 = fullfile('SupplementaryTables', 'Table_S6.xlsx');
SepiAllProteinInfo=struct2table(SepiAllProteinInfo);
CacnesAllProteinInfo=struct2table(CacnesAllProteinInfo);
SepiAnnotationsHomologues=struct2table(SepiAnnotationsHomologues);
CacnesAnnotationsHomologues=struct2table(CacnesAnnotationsHomologues);
cacnes_annotations_all=struct2table(cacnes_annotations_all);
sepi_annotations_all=struct2table(sepi_annotations_all);
cacnes_annotations_all = parse_annotations(cacnes_annotations_all);
sepi_annotations_all = parse_annotations(sepi_annotations_all);
% write to S6
writetable(cacnes_annotations_all, fileS6, 'Sheet', 'cacnes_annotations_all');
writetable(sepi_annotations_all, fileS6, 'Sheet', 'semi_annotations_all');
writetable(HomologousClustersCacnesAssemblies, fileS6, 'Sheet', 'HomologMembershipCacnes');
writetable(HomologousClustersSepiAssemblies, fileS6, 'Sheet', 'HomologMembershipSepi');
writetable(CacnesAllProteinInfo, fileS6, 'Sheet', 'CacnesAllProteinInfo');
writetable(SepiAllProteinInfo, fileS6, 'Sheet', 'SepiAllProteinInfo');
writetable(CacnesAnnotationsHomologues, fileS6, 'Sheet', 'CacnesAnnotationsHomologs');
writetable(SepiAnnotationsHomologues, fileS6, 'Sheet', 'SepiAnnotationsHomologs');
end

%% Local helper function example
function Mnew = parse_to_csv_version(M)
% Convert fields in M to a table format suitable for CSV/Excel output.
Mnew = table;
Mnew.SID = M.SID;
Mnew.TP = M.TP;
if ~isstring(M.location)
    Mnew.location = arrayfun(@(x) strjoin(x{:}','\t'), M.location);
end
if ~isstring(M.TP)
    try
        Mnew.time_points = arrayfun(@(x) strjoin(string(x(:)'),'\t'), M.TP);
    catch
        Mnew.time_points = arrayfun(@(x) strjoin(string(x{:}'),'\t'), M.TP);
    end
    Mnew.TP=[];
end
Mnew.Age = M.Age;
Mnew.Sex = M.Sex;
Mnew.CacnesPhylotypeNames = arrayfun(@(x) strjoin(M.CacnesPhylotypeNames(x,:), '\t'), 1:size(M,1))';
Mnew.SepiPhylotypeNames = repmat('1 2 3 4', size(M,1), 1);
Mnew.SepiSubPhylotypeNames = arrayfun(@(x) strjoin(M.SepiSubphylotypeNames(x,:), '\t'), 1:size(M,1))';
Mnew.SepiLineageNumbers = repmat(strjoin(string(M.SepiLineageNumbers(1,:)), '\t'), size(M,1), 1);
Mnew.CacnesLineageNumbers = repmat(strjoin(string(M.CacnesLineageNumbers(1,:)), '\t'), size(M,1), 1);
Mnew.SpeciesNamesBracken = repmat(strjoin(string(M.NamesSpecies(1,:)), '\t'), size(M,1), 1);
Mnew.GenusNamesBracken = repmat(strjoin(string(M.NamesGenus(1,:)), '\t'), size(M,1), 1);
Mnew.subject_cutotype = M.subject_cutotype;
if any(ismember(M.Properties.VariableNames, "cacnes_subphylolev_concatenated"))
    M.CombinedCacnesSubphylo = M.cacnes_subphylolev_concatenated;
end
R = size(M,1);
Variables = ["SepiMGCoverage", "SepiOverHuman", "CacnesMGCoverage", "CacnesOverHuman", ...
    "FilteredSpeciesLevelReads", "BrackenAbundanceClustering", "BrackenNamesClustering", ...
    "BrackenGenus", "BrackenSpecies", "CacnesLineageAbundance", "CacnesPhylotypeAbundance", ...
    "SepiLineageAbundance", "SepiPhylotypeAbundance", "subject_cutotype", ...
    "umap_clusters_subjecttime_bracken_species", "CombinedCacnesLineages", "CombinedCacnesPhylotypes", ...
    "CombinedSepiLineages", "CombinedSepiSubphylo", "CombinedSepiPhylo", "CombinedCacnesSubphylo"];
for i = 1:numel(Variables)
    try
        col = M.(Variables(i));
    catch
        continue;
    end
    if iscell(col)
        col = vertcat(col{:});
    end
    col_new = arrayfun(@(x) {strjoin(string(col(x,:)), '\t')}, 1:R);
    Mnew.(Variables(i)) = vertcat(col_new{:});
end
end

%% parse lineage dmrca data
function LineageDMRCAsForSupplement = LineageData2supplement(LineageDMRCAsForSupplement)

LineageDMRCAsForSupplement = removevars(LineageDMRCAsForSupplement, ["RTT_clade","RTT_clade_ugenotypes","TreeSNPs","Subjects","UniqueSNPs_genotype"]);
LineageDMRCAsForSupplement = removevars(LineageDMRCAsForSupplement, ["SubjectTreeSNPs_unique_genotypes","SubjectTreeSNPs_isolates"]);

LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates);
LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates=vertcat(LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates{:});
%
LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes);
LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes=vertcat(LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes{:});
%
LineageDMRCAsForSupplement.RTT_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.RTT_subject_isolates);
LineageDMRCAsForSupplement.RTT_subject_isolates=vertcat(LineageDMRCAsForSupplement.RTT_subject_isolates{:});
%
LineageDMRCAsForSupplement.dMRCA_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.dMRCA_subject_isolates);
LineageDMRCAsForSupplement.dMRCA_subject_isolates=vertcat(LineageDMRCAsForSupplement.dMRCA_subject_isolates{:});
%
LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes);
LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes=vertcat(LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes{:});

LineageDMRCAsForSupplement.RTT_subject_ugenotypes=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.RTT_subject_ugenotypes);
LineageDMRCAsForSupplement.RTT_subject_ugenotypes=vertcat(LineageDMRCAsForSupplement.RTT_subject_ugenotypes{:});
end

%% parse gain loss
function anno_out = parse_gain_loss(anno_in)
R = size(anno_in,1);
Variables = string(anno_in.Properties.VariableNames);
V = numel(Variables);

anno_out=table;
for i = 1:V
    col=anno_in.(Variables(i));
    if isempty(col)
        anno_out.(Variables(i))=repmat("",R,1);
        continue
    end
    if iscell(col)
        col_new = arrayfun(@(x) strjoin(string(col{x}),'\t'), 1:R);
    else
        col_new = arrayfun(@(x) strjoin(string(col(x,:)),'\t'), 1:R);
    end
    if mean(col_new=="")==1
        continue
    else
        anno_out.(Variables(i))=col_new';
    end
end
end

%% parse annotations
function anno_out = parse_annotations(anno_in)


R = size(anno_in,1);
Variables = string(anno_in.Properties.VariableNames);
V = numel(Variables);

anno_out=table;
for i = 1:numel(Variables)
    col=anno_in.(Variables(i));
    if iscell(col)
        col_new = arrayfun(@(x) strjoin(string(col{x}),'\t'), 1:R);
    else
        col_new = arrayfun(@(x) strjoin(string(col(x,:)),'\t'), 1:R);
    end
    if mean(col_new=="")==1
        continue
    else
        anno_out.(Variables(i))=col_new';
    end
end

end