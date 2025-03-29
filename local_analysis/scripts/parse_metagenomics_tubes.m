function [MSubject, MsubjectTime, Mtubes]=parse_metagenomics_tubes(M,ConcatenatedReadsData)

% This function takes the table M_all_tubes and parses it into a table for
% each subject, subjecttimepoint, and tube. For tables involving more than
% one tube, it normalizes/averages where approproate and also adds data
% from concatenatedreadsdata

R = size(M,1);
MFields = string(M.Properties.VariableNames);
IndividualTubes = arrayfun(@(x) strjoin([M.SID(x) string(M.TP(x)) M.location(x)]) , 1:R);

[uIndividualTubes,ia,ic]=unique(IndividualTubes);

% these fields will not change for individual tubes
InvariantFields = ["SID" "TP" "location" "Age" "Sex" "CacnesPhylotypeNames" "CacnesLineageNumbers" "CacnesLineagePhylotypes" "SepiLineageNumbers" "SepiLineagePhylotypes" "SepiSubphylotypeNames" "NamesSpecies" "NamesGenus"];
% these fields need to be combined
FieldsToCombine =["tube" "plate" "well" "index" "CacnesMGCoverage" "SepiMGCoverage" "SampleNames" "BrackenGenus" "BrackenReadsAssignedGenus" "BrackenReadsAssignedSpecies" "BrackenSpecies"];

is_InvariantField = arrayfun(@(x) find(x==MFields), InvariantFields);
Mtubes=M(ia,is_InvariantField);

% now add the index for each tube
Mtubes.index604 = arrayfun(@(x) {M.index(IndividualTubes==x)}, uIndividualTubes)';
index = arrayfun(@(x) {find(IndividualTubes==x)}, uIndividualTubes)';

% for each of the variable fields, add it to Mtubes

for i = 1:numel(FieldsToCombine)
    F = FieldsToCombine(i);
    Mtubes.(F) = arrayfun(@(x) {M.(F)(x{:},:)},  index);
end

% for the following fields, add together
Mtubes.CacnesMGCoverage = cellfun(@sum, Mtubes.CacnesMGCoverage);
Mtubes.SepiMGCoverage = cellfun(@sum, Mtubes.SepiMGCoverage);
Mtubes.BrackenReadsAssignedGenus=arrayfun(@(x) {sum(x{:},1)}, Mtubes.BrackenReadsAssignedGenus);
Mtubes.BrackenReadsAssignedSpecies=arrayfun(@(x) {sum(x{:},1)}, Mtubes.BrackenReadsAssignedSpecies);

% For these ones, normalize together
Mtubes.BrackenGenus=cellfun(@(x) {mean(x,1)./sum(mean(x,1))}, Mtubes.BrackenGenus);
Mtubes.BrackenSpecies=cellfun(@(x) {mean(x,1)./sum(mean(x,1))}, Mtubes.BrackenSpecies);

% for metagenomics assignments, pick the one with the most assigned 
%% First, reduce M to unique tubes (some tubes were sequenced twice)
% these data are not currently used in any main figures
Mtubes=add_phlame_data_tubes(Mtubes,M);

%% per subject/timepoint (Combining multiple tubes from a single timepoint)

R = size(Mtubes,1);

SubjectTime = arrayfun(@(x) strjoin([Mtubes.SID(x) string(Mtubes.TP(x))],'_') , 1:R);
MFields = string(Mtubes.Properties.VariableNames);

[uSubjectTime,ia,ic]=unique(SubjectTime);

% these fields will not change for individual tubes
InvariantFields = ["SID" "TP" "Age" "Sex" "CacnesPhylotypeNames" "CacnesLineageNumbers" "CacnesLineagePhylotypes" "SepiLineageNumbers" "SepiLineagePhylotypes" "SepiSubphylotypeNames" "NamesSpecies" "NamesGenus"];
% these fields need to be combined
FieldsToCombine =["location" "tube" "plate" "well" "index" "CacnesMGCoverage" "SepiMGCoverage" "SampleNames" "BrackenGenus" "BrackenReadsAssignedGenus" "BrackenReadsAssignedSpecies" "BrackenSpecies"];

is_InvariantField = arrayfun(@(x) find(x==MFields), InvariantFields);
MsubjectTime=Mtubes(ia,is_InvariantField);

% now add the index for each tube
index = arrayfun(@(x) {find(SubjectTime==x)}, uSubjectTime)';

% for each of the variable fields, add it to MsubjectTime

for i = 1:numel(FieldsToCombine)
    F = FieldsToCombine(i);
    MsubjectTime.(F) = arrayfun(@(x) {Mtubes.(F)(x{:},:)},  index);
end

% for the following fields, add together
MsubjectTime.CacnesMGCoverage = cellfun(@sum, MsubjectTime.CacnesMGCoverage);
MsubjectTime.SepiMGCoverage = cellfun(@sum, MsubjectTime.SepiMGCoverage);
MsubjectTime.BrackenReadsAssignedGenus=arrayfun(@(x) {sum(vertcat(x{:}{:}),1)}, MsubjectTime.BrackenReadsAssignedGenus);
MsubjectTime.BrackenReadsAssignedSpecies=arrayfun(@(x) {sum(vertcat(x{:}{:}),1)}, MsubjectTime.BrackenReadsAssignedSpecies);

% For these ones, normalize together
MsubjectTime.BrackenGenus=cellfun(@(x) {mean(vertcat(x{:}),1)./sum(mean(vertcat(x{:}),1))}, MsubjectTime.BrackenGenus);
MsubjectTime.BrackenSpecies=cellfun(@(x) {mean(vertcat(x{:}),1)./sum(mean(vertcat(x{:}),1))}, MsubjectTime.BrackenSpecies);

MsubjectTime.BrackenSpecies=vertcat(MsubjectTime.BrackenSpecies{:});
MsubjectTime.BrackenGenus=vertcat(MsubjectTime.BrackenGenus{:});
MsubjectTime.BrackenReadsAssignedGenus=vertcat(MsubjectTime.BrackenReadsAssignedGenus{:});
MsubjectTime.BrackenReadsAssignedSpecies=vertcat(MsubjectTime.BrackenReadsAssignedSpecies{:});

% add concatenated reads data
SubjectTime = ~ismissing(ConcatenatedReadsData.TP)&ismissing(ConcatenatedReadsData.Plate);
ConcatenatedReadsDataSubjectTime=ConcatenatedReadsData(SubjectTime,:);

% check to make sure the order is right
idx = arrayfun(@(x) find(x==ConcatenatedReadsDataSubjectTime.FullNames), uSubjectTime);
% Sepi clusters
MsubjectTime.sepi_clusterlev_concatenated_conservative=ConcatenatedReadsDataSubjectTime.sepi_clusterlev_concatenated_conservative(idx,:);
MsubjectTime.sepi_clusterlev_concatenated_loose=ConcatenatedReadsDataSubjectTime.sepi_clusterlev_concatenated_loose(idx,:);
% Cacnes clusters
MsubjectTime.cacnes_clusterlev_concatenated_conservative=ConcatenatedReadsDataSubjectTime.cacnes_clusterlev_concatenated_conservative(idx,:);
MsubjectTime.cacnes_clusterlev_concatenated_loose=ConcatenatedReadsDataSubjectTime.cacnes_clusterlev_concatenated_loose(idx,:);
% Sepi phylotypes
MsubjectTime.sepi_phylolev_concatenated_conservative=ConcatenatedReadsDataSubjectTime.sepi_phylolev_concatenated_conservative(idx,:);
MsubjectTime.sepi_phylolev_concatenated_loose=ConcatenatedReadsDataSubjectTime.sepi_phylolev_concatenated_loose(idx,:);
% cacnes phylotypes
MsubjectTime.cacnes_phylolev_concatenated_loose=ConcatenatedReadsDataSubjectTime.cacnes_phylolev_concatenated_loose(idx,:);
MsubjectTime.cacnes_phylolev_concatenated_conservative=ConcatenatedReadsDataSubjectTime.cacnes_phylolev_concatenated_conservative(idx,:);
% Sepi subphylotypes
MsubjectTime.sepi_subphylolev_concatenated_loose=ConcatenatedReadsDataSubjectTime.sepi_subphylolev_concatenated_loose(idx,:);
MsubjectTime.sepi_subphylolev_concatenated_conservative=ConcatenatedReadsDataSubjectTime.sepi_subphylolev_concatenated_conservative(idx,:);
% Cacnes subphylotypes
MsubjectTime.cacnes_subphylolev_concatenated=ConcatenatedReadsDataSubjectTime.cacnes_subphylolev_concatenated(idx,:);



%% by subject



Subject = Mtubes.SID;
MFields = string(Mtubes.Properties.VariableNames);

[uSubject,ia,ic]=unique(Subject);

% these fields will not change for individual tubes
InvariantFields = ["SID" "Sex" "CacnesPhylotypeNames" "CacnesLineageNumbers" "CacnesLineagePhylotypes" "SepiLineageNumbers" "SepiLineagePhylotypes" "SepiSubphylotypeNames" "NamesSpecies" "NamesGenus"];
% these fields need to be combined
FieldsToCombine =["location" ,"TP", "tube" "Age" "plate" "well" "index" "CacnesMGCoverage" "SepiMGCoverage" "SampleNames" "BrackenGenus" "BrackenReadsAssignedGenus" "BrackenReadsAssignedSpecies" "BrackenSpecies"];

is_InvariantField = arrayfun(@(x) find(x==MFields), InvariantFields);
MSubject=Mtubes(ia,is_InvariantField);

% now add the index for each tube
index = arrayfun(@(x) {find(Subject==x)}, uSubject)';

% for each of the variable fields, add it to MSubject

for i = 1:numel(FieldsToCombine)
    F = FieldsToCombine(i);
    MSubject.(F) = arrayfun(@(x) {Mtubes.(F)(x{:},:)},  index)';
end

% for the following fields, add together
MSubject.CacnesMGCoverage = cellfun(@sum, MSubject.CacnesMGCoverage);
MSubject.SepiMGCoverage = cellfun(@sum, MSubject.SepiMGCoverage);
MSubject.BrackenReadsAssignedGenus=arrayfun(@(x) {sum(vertcat(x{:}{:}),1)}, MSubject.BrackenReadsAssignedGenus);
MSubject.BrackenReadsAssignedSpecies=arrayfun(@(x) {sum(vertcat(x{:}{:}),1)}, MSubject.BrackenReadsAssignedSpecies);
%
MSubject.Age = cellfun(@(x) mean(unique(x)), MSubject.Age);
% For these ones, normalize together
MSubject.BrackenGenus=cellfun(@(x) {mean(vertcat(x{:}),1)./sum(mean(vertcat(x{:}),1))}, MSubject.BrackenGenus);
MSubject.BrackenSpecies=cellfun(@(x) {mean(vertcat(x{:}),1)./sum(mean(vertcat(x{:}),1))}, MSubject.BrackenSpecies);

% add concatenated reads data
Subject = ismissing(ConcatenatedReadsData.TP)&ismissing(ConcatenatedReadsData.Plate);
ConcatenatedReadsDataSubject=ConcatenatedReadsData(Subject,:);

% check to make sure the order is right
idx = arrayfun(@(x) find(x==ConcatenatedReadsDataSubject.FullNames), uSubject);
% Sepi clusters
MSubject.sepi_clusterlev_concatenated_conservative=ConcatenatedReadsDataSubject.sepi_clusterlev_concatenated_conservative(idx,:);
MSubject.sepi_clusterlev_concatenated_loose=ConcatenatedReadsDataSubject.sepi_clusterlev_concatenated_loose(idx,:);
% Cacnes clusters
MSubject.cacnes_clusterlev_concatenated_conservative=ConcatenatedReadsDataSubject.cacnes_clusterlev_concatenated_conservative(idx,:);
MSubject.cacnes_clusterlev_concatenated_loose=ConcatenatedReadsDataSubject.cacnes_clusterlev_concatenated_loose(idx,:);
% Sepi phylotypes
MSubject.sepi_phylolev_concatenated_conservative=ConcatenatedReadsDataSubject.sepi_phylolev_concatenated_conservative(idx,:);
MSubject.sepi_phylolev_concatenated_loose=ConcatenatedReadsDataSubject.sepi_phylolev_concatenated_loose(idx,:);
% cacnes phylotypes
MSubject.cacnes_phylolev_concatenated_loose=ConcatenatedReadsDataSubject.cacnes_phylolev_concatenated_loose(idx,:);
MSubject.cacnes_phylolev_concatenated_conservative=ConcatenatedReadsDataSubject.cacnes_phylolev_concatenated_conservative(idx,:);
% Sepi subphylotypes
MSubject.sepi_subphylolev_concatenated_loose=ConcatenatedReadsDataSubject.sepi_subphylolev_concatenated_loose(idx,:);
MSubject.sepi_subphylolev_concatenated_conservative=ConcatenatedReadsDataSubject.sepi_subphylolev_concatenated_conservative(idx,:);
%Cacnes subphylotypes
MSubject.cacnes_subphylolev_concatenated=ConcatenatedReadsDataSubject.cacnes_subphylolev_concatenated(idx,:);

Mtubes.BrackenGenus=vertcat(Mtubes.BrackenGenus{:});
MSubject.BrackenGenus=vertcat(MSubject.BrackenGenus{:});
