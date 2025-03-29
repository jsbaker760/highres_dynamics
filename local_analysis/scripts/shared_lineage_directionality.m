function [f1, f2]= shared_lineage_directionality(L,n_total,total_shared)
%% Get L where there are more than 2 subjects in a lineage

AnyOtherSubjectsThisLineage = arrayfun(@(x) any(L.LineageNumber(1:(size(L,1))~=x)==L.LineageNumber(x)), 1:size(L,1));
% remove lineages containing only one subject with enough isolates
L=L(AnyOtherSubjectsThisLineage,:);

%% Difference between dMRCA subject and dMRCA of the entire lineages

[dDMRCA,i]=sort(L.dDMRCA,'descend');
L =L(i,:);
R = (size(L,1));
C=L.Cutotype;

%% get the number of lineages with sufficient data

NumberLineagesWithSufficientData = numel(unique(L.LineageNumber));

%% Find out which subjects in lineages are sources
% boolean for each row, true where rows are in the same lineages but not
% the same row
OtherSubjectsSameLineage = arrayfun(@(x) {(L.LineageNumber==L.LineageNumber(x)&(1:numel(dDMRCA)~=x)')}, 1:R);
% cutotypes of other subjects in the same lineage (non-self)
CutotypeOtherSubjectSameLineage=cellfun(@(x) {L.Cutotype(x)}, OtherSubjectsSameLineage);
% same size as above, and true where this subject's MRCA is the root i.e.
% is a possible source
SourceOtherSubjectSameLineage=cellfun(@(x) {L.MRCA_equals_root(x)}, OtherSubjectsSameLineage);
% cutotypes for above
CutotypesSources = arrayfun(@(x) {CutotypeOtherSubjectSameLineage{x}(SourceOtherSubjectSameLineage{x}==1)}, 1:R);

%% Delineate whether sources are only parent, only child, both, or neither (shared unknown source)
NoOtherSubjectsRoot = cellfun(@isempty,CutotypesSources);
ParentSource = cellfun(@(x) all(unique(x)==3)&~isempty(x), CutotypesSources);
ChildSource = cellfun(@(x) all(unique(x)<3)&~isempty(x), CutotypesSources);
ComplexSource = cellfun(@(x) any(unique(x)==3)&any(unique(x)<3)&~isempty(x), CutotypesSources);

%% get colors for dots (based on cutotype)
clrs_dots = zeros(numel(C),3);
clrs_dots(C==1,:)=repmat([.337 .705 .913],sum(C==1),1);
clrs_dots(C==2,:)=repmat([0 .4471 .698],sum(C==2),1);
clrs_dots(C==3,:)=repmat([.8253 .3686 .0039],sum(C==3),1);

%% color arrows based on who is source
lrs_arrows = zeros(numel(C),3);
clrs_arrows(ParentSource,:)=repmat([.8253 .3686 .0039],sum(ParentSource),1);
clrs_arrows(ChildSource,:)=repmat([0 0 1],sum(ChildSource),1);
clrs_arrows(ComplexSource,:)=repmat([.85 .85 .85],sum(ComplexSource),1);
clrs_arrows(NoOtherSubjectsRoot,:)=repmat([.85 .85 .85],sum(NoOtherSubjectsRoot),1);

%% find the lineages with clear directionality
NumberLineagesWithClearDirectionality = numel(unique(L.LineageNumber(~ComplexSource&~NoOtherSubjectsRoot&dDMRCA'>0)));

ddMRCA_clear = dDMRCA(~ComplexSource&~NoOtherSubjectsRoot&dDMRCA'>0);
clrs_arrow_clear = clrs_arrows(~ComplexSource&~NoOtherSubjectsRoot&dDMRCA'>0,:);
clrs_dots_clear = clrs_dots(~ComplexSource&~NoOtherSubjectsRoot&dDMRCA'>0,:);

%% plot all the clear instances

R = numel(ddMRCA_clear);
f2=figure;hold on;
arrayfun(@(x) plot([x x],[0 ddMRCA_clear(x)],'LineWidth',1,'Color',clrs_arrow_clear(x,:)) ,1:R)

scatter(1:R,ddMRCA_clear,'filled','CData',clrs_dots_clear)
pbaspect([2,1,1])
xlabel('rank')
ylabel('∆MRCA')
xlim([0 R])

%% make a stacked bar chart

SharedDirectionalClear = NumberLineagesWithClearDirectionality;
SharedUnclear = total_shared-SharedDirectionalClear;
Unshared = n_total-SharedDirectionalClear-SharedUnclear;
f1=figure;
b=bar(1,[SharedDirectionalClear SharedUnclear Unshared],'stacked','FaceColor','flat','BarWidth',.5);
b(1).CData=[0 168 117]./255;
b(2).CData=[83 192 168]./255;
b(3).CData=[.85 .85 .85];

l=legend;
l.String=["has clear directionality" "unclear direction" " unshared"];

xlim([.75 1.25])
ylim([0 90]);
yticks(0:10:90)
pbaspect([1 2 1])
l.Location="westoutside";

[char(string(SharedDirectionalClear)) '/' char(string(NumberLineagesWithSufficientData)) ' lineages with sufficient data have clear directionality']
end