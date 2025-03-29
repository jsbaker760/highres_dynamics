function [NewNumbers] =  getNewClusterIDs(oldclusters,newclusters)
%% initialize
N = max(newclusters);
O = max(oldclusters);
%% relate indices
NewInOld = arrayfun(@(x) {newclusters(oldclusters==x)}, 1:O);
idxold = arrayfun(@(x) {find(oldclusters==x)}, 1:O);
OldInNew    = arrayfun(@(x) {oldclusters(newclusters==x)}, 1:N);
idxnew = arrayfun(@(x) {find(newclusters==x)}, 1:N);

newclusters_inold = arrayfun(@(x) {unique(x{:}(x{:}>0))}, NewInOld);
oldclusters_innew = arrayfun(@(x) {unique(x{:}(x{:}>0))}, OldInNew);

%% find clusters needed to split names
need2split = find(cellfun(@numel,newclusters_inold)>1);
if any(need2split)
    NewInOld_split = NewInOld{need2split};
    idxold_split = idxold{need2split};
    oldclusters(idxold_split(NewInOld_split==109))=73;
    N = max(newclusters);
    O = max(oldclusters);
    NewInOld = arrayfun(@(x) {newclusters(oldclusters==x)}, 1:O);
    idxold = arrayfun(@(x) {find(oldclusters==x)}, 1:O);
    OldInNew    = arrayfun(@(x) {oldclusters(newclusters==x)}, 1:N);
    idxnew = arrayfun(@(x) {find(newclusters==x)}, 1:N);
    newclusters_inold = arrayfun(@(x) {unique(x{:}(x{:}>0))}, NewInOld);
    oldclusters_innew = arrayfun(@(x) {unique(x{:}(x{:}>0))}, OldInNew);
end

%% get new numbers 
[NewNumbers]=deal(zeros(N,1));

if all(cellfun(@numel,newclusters_inold)<=1)&all(cellfun(@numel,oldclusters_innew)<=1) &...
        numel(unique(horzcat(newclusters_inold{:})))==numel((unique(horzcat(newclusters_inold{:})))) & ...
        numel(unique(horzcat(oldclusters_innew{:})))==numel((unique(horzcat(oldclusters_innew{:}))))
    StartIdxNewClusters = O+1;
    for n=1:N
        match = unique(OldInNew{n});
        match = match(match>0);
        if ~isempty(match)
            NewNumbers(idxnew{n})=match;
        else
            NewNumbers(idxnew{n})=StartIdxNewClusters;
            StartIdxNewClusters =StartIdxNewClusters+1;
        end

    end
N = max(NewNumbers);
O = max(newclusters);

NewInOld = arrayfun(@(x) {NewNumbers(newclusters==x)}, 1:O);
idxold = arrayfun(@(x) {find(newclusters==x)}, 1:O);
idxnew = arrayfun(@(x) {find(NewNumbers==x)}, 1:N);
newclusters_inold = arrayfun(@(x) unique(x{:}(x{:}>0)), NewInOld);

for o=1:O
    if ~all(idxnew{newclusters_inold(o)}'==idxold{o})
        error('the idices do not perfectly match')
    end
end
end