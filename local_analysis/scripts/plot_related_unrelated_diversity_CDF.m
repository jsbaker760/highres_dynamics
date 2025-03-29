function [f,f2]=plot_related_unrelated_diversity_CDF(M,cutoff_include,cutoff_real)

f=figure';
%% C. acnes lineage level
subplot(2,1,1)
BoolInclude = sum(M.CombinedCacnesLineages,2)>cutoff_include;
[~]=plot_intra_inter1(M(BoolInclude,:),M.CombinedCacnesLineages(BoolInclude,:),cutoff_real,'C. acnes lineage level');

%% S. epi lineage level
subplot(2,1,2)
BoolInclude = sum(M.CombinedSepiLineages,2)>cutoff_include;
[~]=plot_intra_inter1(M(BoolInclude,:),M.CombinedSepiLineages(BoolInclude,:),cutoff_real,'S. epi lineage level');

%% Supplemental data
% C. acnes lineage level
f2=figure;
subplot(2,1,1)
BoolInclude = sum(M.CombinedCacnesLineages,2)>cutoff_include;
[~]=plot_intra_inter1_supp(M(BoolInclude,:),M.CombinedCacnesLineages(BoolInclude,:),cutoff_real,'C. acnes lineage level');
% S. epi lineage level
subplot(2,1,2)
BoolInclude = sum(M.CombinedSepiLineages,2)>cutoff_include;
[~]=plot_intra_inter1_supp(M(BoolInclude,:),M.CombinedSepiLineages(BoolInclude,:),cutoff_real,'S. epi lineage level');

%% pull out the distributions
function Distributions =plot_intra_inter1(JaccardSubjects,A,cutoffReal,str)
%% get unique subject

SID=JaccardSubjects.SID;
uSID = unique(SID);
SIDidx = arrayfun(@(x) {find(x==SID)}, uSID);

%% get jaccard of each row for self over time comparison

JaccardSubjectTime = jaccard(A);
SameSubjectDifferentTime = (SID==SID')&(JaccardSubjects.TP~=JaccardSubjects.TP');
UniqueEntry = tril(ones(size(JaccardSubjectTime)),-1)==1;
JaccardSameSubjectDifferentTime = JaccardSubjectTime(SameSubjectDifferentTime&UniqueEntry);


%% find every lineage found at any timepoint

% true where lineage is ever found above detection limit
A = arrayfun(@(x) {sum(A(x{:},:)>cutoffReal,1)>0}, SIDidx);
LineageEverFound=vertcat(A{:});

% jaccard distance
JaccardSubjects = jaccard(LineageEverFound);
% get booleans for comparing distances of different groups

FamilyNumber = get_family_numbers(uSID);
NonSelf = tril(ones(size(JaccardSubjects)),-1)==1;
SameFamily = FamilyNumber==FamilyNumber';
SameSubject = uSID==uSID';


%% split distance matrix into different groups to compare



group_names=["Unrelated subjects, any timepoint" "Related subjects, any timepoint" "Same subject, different timepoint"];
Distributions = [{JaccardSubjects(NonSelf&~SameFamily)} {JaccardSubjects(NonSelf&SameFamily&~SameSubject)} {JaccardSameSubjectDifferentTime}];
% turn to 1D arrays for swarmchart
clrs = [{[.85 .85 .85]} {[0 1 1]} {[1 0 1]}];

NonFamilySharing = Distributions{1};
[i,j]=find(NonSelf&~SameFamily);
NonZeroSharing = NonFamilySharing~=1;

['Number of sharing events between non-family members = ' char(string(sum(NonZeroSharing)))]

%%

for i = 3:-1:1
    c=cdfplot(Distributions{i});
    c.Color = clrs{i};
    c.LineWidth=2;
    grid off
    hold on
end

pbaspect([1 1 1])
p12 = ranksum(Distributions{1},Distributions{2});
p13 = ranksum(Distributions{1},Distributions{3});
p23 = ranksum(Distributions{2},Distributions{3});

yticks([0:.25:1])
xticks([0:.25:1])
text(1.1,1.1,['P Unrelated-Related =' char(string(round(p12,2,'significant')))])
text(1.1,1,['P Related-SameSubject =' char(string(round(p23,2,'significant')))])
text(1.1,.9,['P UnrelatedSameSubject =' char(string(round(p13,2,'significant')))])

title(str)
xlabel('Bray-Curtis Dissimilarity')
ylabel('cumulative probability')

function Distributions=plot_intra_inter1_supp(JaccardSubjects,A,cutoffReal,str)
%% get unique subject
Distributions=[];

SID=JaccardSubjects.SID;
uSID = unique(SID);
SIDidx = arrayfun(@(x) {find(x==SID)}, uSID);

%% get jaccard of each row for self over time comparison

JaccardSelf = jaccard(A);
SameSubjectDifferentTime = (SID==SID')&(JaccardSubjects.TP~=JaccardSubjects.TP');
UniqueEntry = tril(ones(size(JaccardSelf)),-1)==1;
Parent = contains(JaccardSubjects.SID,"P");
ParentParent_tp = Parent&Parent';
ChildChild_tp = ~Parent&~Parent';

%% find every lineage found at any timepoint

% true where lineage is ever found above detection limit
LineageEverFound = arrayfun(@(x) {sum(A(x{:},:)>cutoffReal,1)>0}, SIDidx);
LineageEverFound=vertcat(LineageEverFound{:});

% jaccard distance
JaccardSubjects = jaccard(LineageEverFound);
Parent = contains(uSID,"P");
ParentParent_subject = Parent&Parent';
ChildChild_subject = ~Parent&~Parent';
ParentChild_subject= (~Parent&Parent')|(Parent&~Parent');

%% get booleans for comparing distances of different groups
FamilyNumber = get_family_numbers(uSID);
NonSelf = tril(ones(size(JaccardSubjects)),-1)==1;
SameFamily = FamilyNumber==FamilyNumber';

%% split distance matrix into different groups to compare

Distributions_self_time = [{JaccardSelf(SameSubjectDifferentTime&UniqueEntry&ParentParent_tp)} {JaccardSelf(SameSubjectDifferentTime&UniqueEntry&ChildChild_tp)} ];

hold on
c1=cdfplot(Distributions_self_time{1});
c1.LineStyle='--';
c1.LineWidth=1;
c1.Color='magenta';
%
c1=cdfplot(Distributions_self_time{2});
c1.LineStyle='-';
c1.LineWidth=1;
c1.Color='magenta';
grid off
%

Distributions_non_self_time = [{JaccardSubjects(NonSelf&ParentParent_subject&SameFamily)} {JaccardSubjects(NonSelf&ChildChild_subject&SameFamily)} {JaccardSubjects(ParentChild_subject&SameFamily)}];

hold on
c1=cdfplot(Distributions_non_self_time{1});
c1.LineStyle='--';
c1.LineWidth=1;
c1.Color='green';
%
c1=cdfplot(Distributions_non_self_time{2});
c1.LineStyle='-';
c1.LineWidth=1;
c1.Color='green';
grid off
%
%
c1=cdfplot(Distributions_non_self_time{3});
c1.LineStyle='-.';
c1.LineWidth=1;
c1.Color='green';
grid off

title(str)
l=legend;
l.String = ["Self over time (parents)" "Self over time Children)" "Parents vs. Parents" "Children vs. Children" "Parents vs. Children"];
xlabel('Bray-Curtis Dissimilarity')
ylabel('cumulative probability')
