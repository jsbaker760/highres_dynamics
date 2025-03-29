function [SubjectClusterPairs,MRCAs,MRCA2,anc_nts]=sweeps(LineageData,min_niso)

% No molecular clock
% Most subjects are only  going to have 2 tps
% Use subjects with >2tps on a per-isolate basis to get a better p value
% For each lineage
% corr( root to tip distances, time since initial sampling)
% Different lineages (phylotypes) can be different
% Minor subjects (intralineage) can have different rates of evolution
% Exploration on stability of lineages:Remove minor subjects if they're obscuring

%% Parse data

IT = LineageDataToTable(LineageData);

%% Get pairs

SubjectClusterPairs = arrayfun(@(x) strjoin([IT.Subject(x) IT.SpeciesString(x) string(IT.ClusterString(x))])   , 1:size(IT,1));
[SubjectClusterPairs ,~ ,uidx] = unique(SubjectClusterPairs);

%% run through and pull out info

C = numel(SubjectClusterPairs);
[anc_nts,MRCAs,MRCA2] = deal(cell(1,C));
i = 0;
for c= 1:C
    SubjectTimesAll = (IT.SubjectTime(uidx==c));
    Subject= unique((IT.Subject(uidx==c)));
    uSubjectTimesAll=unique(SubjectTimesAll);
    if numel(uSubjectTimesAll)<2
        continue
    end
    T = tabulate(SubjectTimesAll);
    if ~(T{end,2}>=min_niso&sum([T{1:(end-1),2}])>=min_niso)
        continue
    end
    ldx = find(LineageData.cladenumber==unique(IT.ClusterString(uidx==c))&LineageData.SpeciesName==unique(IT.SpeciesString(uidx==c)));
    anc_nts{c}=LineageData.anc_nti{ldx}(LineageData.goodpos{ldx});
    mutations_all = LineageData.tree_snps{ldx};%(:,~LineageData.outgroup{ldx}==1);
    samplenames=LineageData.samplenames{ldx}(~LineageData.outgroup{ldx}==1);
    mutations_subject=mutations_all(:,startsWith(samplenames,Subject));
    mutations_subject(mean(mutations_subject,2)==1,:)=[];
    samplenames_subject=samplenames(startsWith(samplenames,Subject));   
    is_latest_tp = startsWith(samplenames_subject,uSubjectTimesAll(end));
    muts_earlier = mutations_subject(:,~is_latest_tp);
    muts_last = mutations_subject(:,is_latest_tp);

    MRCAs{c}=muts_earlier;
    MRCA2{c}=muts_last;
end


