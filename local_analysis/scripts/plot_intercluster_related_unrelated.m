function plot_intercluster_related_unrelated(IT,CacnesLineageNumbers,SepiLineageNumbers,CacnesDM,SepiDM,CacnesLineages,SepiLineages)

FamiliesCacnes = arrayfun(@(x) { unique(IT.Family(IT.SpeciesString=="cacnes"&IT.ClusterString==x)) }, CacnesLineageNumbers);
FamiliesSepi = arrayfun(@(x) { unique(IT.Family(IT.SpeciesString=="sepi"&IT.ClusterString==x)) }, SepiLineageNumbers);

AnySameSubjectBoolCacnes = zeros(87);
AnySameSubjectBoolSepi = zeros(76);

for i = 1:87
    for j = 1:87
        AnySameSubjectBoolCacnes(i,j)=numel(intersect(FamiliesCacnes{i},FamiliesCacnes{j}))>0;
    end
end

for i = 1:76
    for j = 1:76
        AnySameSubjectBoolSepi(i,j)=numel(intersect(FamiliesSepi{i},FamiliesSepi{j}))>0;
    end
end

CacnesClusterDM=intercluster_dm_pairs(CacnesDM,CacnesLineages.clusters,CacnesLineageNumbers);
SepiClusterDM=intercluster_dm_pairs(SepiDM,SepiLineages.clusters,SepiLineageNumbers);

f=figure;
subplot(4,1,1)
y = CacnesClusterDM;
b = AnySameSubjectBoolCacnes==1;
bool = tril(ones(size(b)),-1)==1;
Y = y(bool);
is_same = b(bool);
[Ysorted,idx]=sort(Y,'ascend');
is_same_sorted = is_same(idx);
%
x = 1:numel(Ysorted);
bar(x(is_same_sorted),Ysorted(is_same_sorted),'FaceColor','red','EdgeColor','none','BarWidth',1)
hold on
bar(x(~is_same_sorted),Ysorted(~is_same_sorted),'FaceColor','black','EdgeColor','none','BarWidth',1)
xlim([0 100])
xlabel('rank (ascending order)')
ylabel('pairwise inter-lineage distance')
l=legend;
l.String=[{'same family'} {'different family'}];
title('C. acnes inter-cluster distances')

subplot(4,1,2)
bar(x(is_same_sorted),Ysorted(is_same_sorted),'FaceColor','red','EdgeColor','none','BarWidth',1)
hold on
bar(x(~is_same_sorted),Ysorted(~is_same_sorted),'FaceColor','black','EdgeColor','none','BarWidth',1)
xlabel('rank (ascending order)')
ylabel('pairwise inter-lineage distance')

%
l=legend;
l.String=[{'same family'} {'different family'}];
%;
subplot(4,1,3)
y = SepiClusterDM;
b = AnySameSubjectBoolSepi==1;
bool = tril(ones(size(b)),-1)==1;
Y = y(bool);
is_same = b(bool);
[Ysorted,idx]=sort(Y,'ascend');
is_same_sorted = is_same(idx);
%
x = 1:numel(Ysorted);
bar(x(is_same_sorted),Ysorted(is_same_sorted),'FaceColor','red','EdgeColor','none','BarWidth',1)
hold on
bar(x(~is_same_sorted),Ysorted(~is_same_sorted),'FaceColor','black','EdgeColor','none','BarWidth',1)
xlim([0 100])
xlabel('rank (ascending order)')
ylabel('pairwise inter-lineage distance')
l=legend;
l.String=[{'same family'} {'different family'}];

title('S. epidermidis inter-cluster distances')
subplot(4,1,4)
bar(x(is_same_sorted),Ysorted(is_same_sorted),'FaceColor','red','EdgeColor','none','BarWidth',1)
hold on
bar(x(~is_same_sorted),Ysorted(~is_same_sorted),'FaceColor','black','EdgeColor','none','BarWidth',1)
xlabel('rank (ascending order)')
ylabel('pairwise inter-lineage distance')

%
l=legend;
l.String=[{'same family'} {'different family'}];