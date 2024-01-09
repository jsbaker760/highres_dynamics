function f=plot_related_unrelated_diversity(M,cutoff_include,cutoff_real)



f=figure

%% C. acnes lineage level
subplot(2,2,1)
BoolInclude = sum(M.CombinedCacnesLineages,2)>cutoff_include;
plot_intra_inter(M(BoolInclude,:),M.CombinedCacnesLineages(BoolInclude,:),cutoff_real,'C. acnes lineage level')

%% C. acnes phylotype level
subplot(2,2,2)
BoolInclude = sum(M.CombinedCacnesPhylotypes,2)>cutoff_include;
plot_intra_inter(M(BoolInclude,:),M.CombinedCacnesPhylotypes(BoolInclude,:),cutoff_real,'C. acnes phylotype level')

%% S. epi lineage level
subplot(2,2,3)
BoolInclude = sum(M.CombinedSepiLineages,2)>cutoff_include;
plot_intra_inter(M(BoolInclude,:),M.CombinedSepiLineages(BoolInclude,:),cutoff_real,'S. epi lineage level')
%% S. epi phylotype level

subplot(2,2,4)
BoolInclude = sum(M.CombinedSepiPhylo,2)>cutoff_include;
plot_intra_inter(M(BoolInclude,:),M.sepi_phylolev_concatenated_conservative(BoolInclude,:),cutoff_real,'S. epi phylo level')

%%



function p =plot_intra_inter(JaccardSubjects,A,cutoffReal,str)
%% get unique subject

M=JaccardSubjects;

SID=JaccardSubjects.SID;
[uSID,ia,ic] = unique(SID);
SIDidx = arrayfun(@(x) {find(x==SID)}, uSID)



%% get jaccard of each row for self over time comparison



JaccardSelf = jaccard(A);
SameSubjectDifferentTime = (SID==SID')&(JaccardSubjects.TP~=JaccardSubjects.TP');
JaccardSameSubjectDifferentTime = JaccardSelf(SameSubjectDifferentTime);


%% find every lineage found at any timepoint



% true where lineage is ever found above detection limit
LineageEverFound = arrayfun(@(x) {sum(A(x{:},:)>cutoffReal,1)>0}, SIDidx);
LineageEverFound=vertcat(LineageEverFound{:});

% jaccard distance
JaccardSubjects = jaccard(LineageEverFound);
%%

FamilyNumbers = get_family_numbers(uSID)
family_with = LineageEverFound.*FamilyNumbers'
Nfam_with=arrayfun(@(x) numel(unique(family_with(:,x)))>2, 1:size(family_with,2))
idx_with = find(Nfam_with);


%% get booleans for comparing distances of different groups



FamilyNumber = get_family_numbers(uSID);
IsParent = contains(uSID,"P");
NonSelf = tril(ones(size(JaccardSubjects)),-1)==1;
SameFamily = FamilyNumber==FamilyNumber';
ParentParent = IsParent&IsParent';
ParentChild = (IsParent&~IsParent')|(~IsParent&IsParent');
ChildChild = (~IsParent&~IsParent');



%% split distance matrix into different groups to compare



group_names=["Parent-Child(related)" "Parent-Child(unrelated)" "Child-Child(related)" "Child-Child(unrelated)" "Parent-Parent(related)" "Parent-Parent(unrelated)"];
% corrosponds to group names above
Ys=[{JaccardSubjects(NonSelf&SameFamily&ParentChild)} {JaccardSubjects(NonSelf&~SameFamily&ParentChild)} {JaccardSubjects(NonSelf&SameFamily&ChildChild)} {JaccardSubjects(NonSelf&~SameFamily&ChildChild)} {JaccardSubjects(NonSelf&SameFamily&ParentParent)} {JaccardSubjects(NonSelf&~SameFamily&ParentParent)}];
Xs = arrayfun(@(x) {x*ones(size(Ys{x}))}, 1:6);

% turn to 1D arrays for swarmchart
Xs=vertcat(Xs{:});
Ys=vertcat(Ys{:});

% get colormap-- gray is unrelated, cyan is related
related = mod(Xs,2)==1;
clrs = .85*ones(numel(Xs),3);
clrs(related,:)=repmat([0 1 1],sum(related),1);

%%


uFamilies = unique(FamilyNumber);
[TotalLineages,SharedLineages]=deal(zeros(numel(uFamilies),1));
for i = 1:numel(uFamilies)
    SidsFamily = uSID(FamilyNumber==uFamilies(i))
    if sum(contains(SidsFamily,"P"))==2
        ParentSids = SidsFamily(contains(SidsFamily,"P"));
        Lineages = LineageEverFound(ismember(uSID,ParentSids),:);
        TotalLineages(i)=sum(sum(Lineages,1)>0);
        SharedLineages(i)=sum(sum(Lineages,1)==2);
    end
end
TotalLineages = sum(TotalLineages)
SharedLineages=sum(SharedLineages)

proportionshared = SharedLineages/TotalLineages


%% add self data



XsSelf = repmat(7,1,numel(JaccardSameSubjectDifferentTime));
clrs_self = repmat([1 0 1],numel(JaccardSameSubjectDifferentTime),1);

Xs=[Xs' XsSelf];
Ys=[Ys' JaccardSameSubjectDifferentTime'];
clrs=[clrs ; clrs_self];
group_names=[group_names "Same Subject; Different TP"];

%%



%% make swarmchart and label



swarmchart(Xs,Ys,'filled','XJitterWidth',.5,'XJitter','density','CData',clrs)
xticks(1:7)
xticklabels(group_names)
ylim([0 1])
xlim([.5 7.5])
ylabel('Jaccard')
pbaspect([1 .5 1])
title(str)

% add p-values
p12 = ranksum(Ys(Xs==1),Ys(Xs==2));
p34 = ranksum(Ys(Xs==3),Ys(Xs==4));
p56 = ranksum(Ys(Xs==5),Ys(Xs==6));

text(1,1,string(p12))
text(3,1,string(p34))
text(5,1,string(p56))

% add median lines
hold on 
arrayfun(@(x) plot([x-.25 x+.25],[median(Ys(Xs==x)) median(Ys(Xs==x))]), 1:7)