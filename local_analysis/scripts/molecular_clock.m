function [SubjectClusterPairs,tps,root_to_tip]=molecular_clock(LineageData,min_niso)

%% Parse lineage-level data
em = cellfun(@isempty, LineageData.anc_nti);
LineageData(em,:)=[];
IT = LineageDataToTable(LineageData);

%% Get pairs of subject in lineages

SubjectClusterPairs = arrayfun(@(x) strjoin([IT.Subject(x) IT.SpeciesString(x) string(IT.ClusterString(x))])   , 1:size(IT,1));
[SubjectClusterPairs ,~ ,uidx] = unique(SubjectClusterPairs);

%% run through and pull out info

C = numel(SubjectClusterPairs);
[root_to_tip,tps,mRCAdiff] = deal(cell(1,C));
rtt_all = cell(1,C);
i = 0;
for c= 1:C
    SubjectTimesAll = (IT.SubjectTime(uidx==c));
    Subject= unique((IT.Subject(uidx==c)));
    uSubjectTimesAll=unique(SubjectTimesAll);
    if numel(uSubjectTimesAll)<2
        continue
    end
    T = tabulate(SubjectTimesAll);
    good = vertcat(T{:,2})>=min_niso;
    if sum(good)<2
        continue
    end
    uSubjectTimesAll = uSubjectTimesAll(good);
    Tps = arrayfun(@(x) str2double(x{:}(end)), cellstr(unique(uSubjectTimesAll)));
    T = numel(Tps);
    instance_root_to_tip = cell(T,1);
    tps{c}=Tps;
    ldx = find(LineageData.cladenumber==unique(IT.ClusterString(uidx==c))&LineageData.SpeciesName==unique(IT.SpeciesString(uidx==c)));
    mutations_all = LineageData.tree_snps{ldx}(:,~LineageData.outgroup{ldx}==1);
    mutations_all(sum(mutations_all,2)==0|sum(mutations_all,2)==size(mutations_all,2),:)=[];
    samplenames=LineageData.samplenames{ldx}(~LineageData.outgroup{ldx}==1);
    mutations_subject=mutations_all(:,startsWith(samplenames,Subject));
    rtt_all{c}=sum(mutations_subject,1);
    samplenames_subject=samplenames(contains(samplenames,Subject));
    diff = sum(sum(mutations_subject,2)==size(mutations_subject,2));
    mRCAdiff{c}=repmat(diff,T,1);
    for t = 1:T
        subjecttimebool = startsWith(samplenames_subject,uSubjectTimesAll(t));
        instance_root_to_tip{t} = sum(mutations_subject(:,subjecttimebool));
        
    end
    root_to_tip{c}=instance_root_to_tip;
end


