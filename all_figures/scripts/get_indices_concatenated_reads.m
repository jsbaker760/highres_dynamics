function T=get_indices_concatenated_reads(M)
%%

Subject = M.SID;
SubjectTime      = arrayfun(@(x) strjoin([M.SID(x) "TP" string(M.TP(x)) "all"],'_'), 1:size(M,1));
SubjectTimePlate      = arrayfun(@(x) strjoin([M.SID(x) string(M.TP(x)) M.plate(x)],'_'), 1:size(M,1)) ;

%%

FaceSite = ismember(M.location,["F","N","K","C"]);
T = table;
idx = 1;

% first, all the ones by subject
UniqueSubjects = unique(Subject);

%%

for s = 1:numel(UniqueSubjects)
    ind = (M.index(FaceSite&Subject==UniqueSubjects(s)))';
    if ~isempty(ind)
        T.IndicesONEIndexed(idx) = {ind};
        T.SampleName(idx) = strjoin(["ConcatenatedReads" UniqueSubjects(s)],'_');
        idx = idx+1;
    else 
        foo = 1;
    end
    
end
% first, all the ones by subject
UniqueSubjectTime = unique(SubjectTime);
for s = 1:numel(UniqueSubjectTime)
    ind=(M.index(FaceSite'&SubjectTime==UniqueSubjectTime(s)))';
    if ~isempty(ind)
        T.SampleName(idx) = strjoin(["ConcatenatedReads" UniqueSubjectTime(s)],'_');
        T.IndicesONEIndexed(idx) = {ind};
        idx=idx+1;
    else 
        foo = 1
    end
end
% first, all the ones by subject
UniqueSubjectTimePlate = unique(SubjectTimePlate);
for s = 1:numel(UniqueSubjectTimePlate)
    ind = (M.index(FaceSite'&SubjectTimePlate==UniqueSubjectTimePlate(s)))';
    if ~isempty(ind)
        T.SampleName(idx) = strjoin(["ConcatenatedReads" UniqueSubjectTimePlate(s)],'_');
        T.IndicesONEIndexed(idx) = {ind};
        idx=idx+1;
    else
        foo = 1
    end
end
%%

%%
indices = arrayfun(@(x) strjoin(string(T.IndicesONEIndexed{x}),' '),1:size(T,1));
T.IndicesONEIndexed =indices';

%%

T.SampleName = string(T.SampleName)

%%


