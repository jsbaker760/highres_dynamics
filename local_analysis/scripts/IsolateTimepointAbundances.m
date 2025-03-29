function [IsolateArrays] = IsolateTimepointAbundances(IT,Mst)

%% get subject_timepoints


SIDTP = IT.SubjectTime;
CacnesLineageNumbers = unique(IT.ClusterString(IT.SpeciesString=="cacnes"));
SepiLineageNumbers = unique(IT.ClusterString(IT.SpeciesString=="sepi"));


%%  MstCutotypes is the type for the subject


SubjectCutotypes=Mst.subject_cutotype;
SubjectCutotypes(Mst.Age>18)=3;


%% parse into arrays


[uSubjectTime,ia,~]=unique(SIDTP);

% number of rows
R =numel(uSubjectTime);

% get lineage numbers from Mst

% initialize abundances
Acacnes = zeros(size(R,1),numel(CacnesLineageNumbers));
ASepi = zeros(size(R,1),numel(SepiLineageNumbers));
% iterate and pull out numbers of isolates
for i = 1:R
    Acacnes(i,:) = arrayfun(@(x) sum(IT.SubjectTime==uSubjectTime(i)&IT.SpeciesString=="cacnes"&IT.ClusterString==x), CacnesLineageNumbers);
    ASepi(i,:) = arrayfun(@(x) sum(IT.SubjectTime==uSubjectTime(i)&IT.SpeciesString=="sepi"&IT.ClusterString==x), SepiLineageNumbers);
end

% make structure
IsolateArrays=table;
IsolateArrays.SubjectTime=uSubjectTime;
IsolateArrays.CacnesIsolates=Acacnes;
IsolateArrays.SepiIsolates=ASepi;

IsolateArrays.CacnesLineageNumbers=repmat(CacnesLineageNumbers',R,1);
IsolateArrays.SepiLineageNumbers=repmat(SepiLineageNumbers',R,1);
%

TPs=char(IT.SubjectTime);
TPs=str2double(string(TPs(:,4)));
IsolateArrays.TP=TPs(ia);
IsolateArrays.SID=IT.Subject(ia);
% the type for the subject, not the timepoint
IsolateArrays.subject_cutotype=arrayfun(@(x) unique(SubjectCutotypes(x==Mst.SID)), IsolateArrays.SID);
% add age where available
IsolateArrays.Age=arrayfun(@(x) max([Mst.Age(Mst.SID==IsolateArrays.SID(x)&Mst.TP==IsolateArrays.TP(x)) 0]), 1:size(IsolateArrays,1))';
% for a few samples, we don't have the exact age handy
% so we add it by looking for a different timepoint with age
% and then adding the time difference

for i = 1:size(IsolateArrays,1)
    if IsolateArrays.Age(i)==0
        ThisTP=IsolateArrays.TP(i);
        idxOtherTP=find(Mst.SID==IsolateArrays.SID(i),1);
        OTherTP=Mst.TP(idxOtherTP);
        AgeOTherTP=Mst.Age(idxOtherTP);
        timediff=tps2duration(ThisTP,OTherTP);
        if OTherTP>ThisTP
            IsolateArrays.Age(i)=AgeOTherTP-timediff;
        else
            IsolateArrays.Age(i)=AgeOTherTP+timediff;
        end
    end
end