function [IsolateTable,LD] = LineageDataToTable(LD)

%% find lineages to exclude
irrelevant_lineage = arrayfun(@(x) any(startsWith(string(LD.samplenames{x}(~LD.outgroup{x})),"0")), 1:size(LD,1));
LD=LD(~irrelevant_lineage,:);


%%  Reshape 
Nclusters = size(LD,1);

SampleNames = arrayfun(@(x) {LD.samplenames{x}(~LD.outgroup{x})}, 1:Nclusters);
ClusterString = arrayfun(@(x) {repmat(LD.cladenumber(x),numel(SampleNames{x}),1)}, 1:Nclusters);
SpeciesString = arrayfun(@(x) {repmat(LD.SpeciesName(x),numel(SampleNames{x}),1)}, 1:Nclusters);

IsolateTable = table;
IsolateTable.SampleNames= string(horzcat(SampleNames{:}))';
IsolateTable.ClusterString= vertcat(ClusterString{:});
IsolateTable.SpeciesString= vertcat(SpeciesString{:});


%% Add metadata
SubjectTimeLocation = arrayfun(@(x) {strsplit((x),'_')}, IsolateTable.SampleNames);

SubjectTimeLocation=arrayfun(@(x) x{:}(1),  SubjectTimeLocation);
Char = cellstr(SubjectTimeLocation);
SubjectTime = arrayfun(@(x) string(x{:}(1:end-1)), SubjectTimeLocation);
Subject = arrayfun(@(x) string(x{:}(1:end-2)), SubjectTimeLocation);
IsolateTable.SubjectTimeLocation=SubjectTimeLocation;  
IsolateTable.SubjectTime=SubjectTime;
IsolateTable.Subject=Subject;
IsolateTable.Family=arrayfun(@(x) {str2double(extract(x,regexpPattern('^[0-9]+')))}, IsolateTable.SubjectTimeLocation);
IsolateTable.Family(cellfun(@numel, IsolateTable.Family)==0)={0};
IsolateTable.Family =vertcat(IsolateTable.Family{:});

end

