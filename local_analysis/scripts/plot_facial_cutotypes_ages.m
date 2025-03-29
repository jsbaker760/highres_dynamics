function plot_facial_cutotypes_ages(MsubjectTime,Msubject,ManuscriptColors)


FC = MsubjectTime.clusters_subjecttime_bracken_species;
Age = MsubjectTime.Age;

f=figure;
clrs2 = [0 0 0; 1 0 0];
subplot (3,1,1)
swarmchart(FC,Age,'filled','CData',clrs2(FC,:))
xlim([0 3])
hold on
arrayfun(@(x) plot([x-.4 x+.4],[median(Age(FC==x)) median(Age(FC==x))],'Color','black','LineWidth',1), 1:2)
pbaspect([1 1 1])
xticks(1:2)
xticklabels(["FC1" "FC2"])
ylabel('Age at sampling (years)')
[p] = ranksum(Age(FC==1),Age(FC==2));
text(1,65,string(p))
%
subplot (3,1,2)
FC = MsubjectTime.ClustersSubject;
Age = MsubjectTime.Age;
[~,idx] = unique(MsubjectTime.SID);
FC =FC(idx);
Age = Age(idx);
swarmchart(FC,Age,'filled','CData',clrs2(FC,:))
xlim([0 3])
hold on
arrayfun(@(x) plot([x-.4 x+.4],[median(Age(FC==x)) median(Age(FC==x))],'Color','black','LineWidth',1), 1:2)
pbaspect([1 1 1])
xticks(1:2)
xticklabels(["FC1" "FC2"])
ylabel('Age at sampling (years)')
[p] = ranksum(Age(FC==1),Age(FC==2));
text(1,65,string(p))
%
subplot (3,1,3)
FC = Msubject.subject_cutotype;
Age = Msubject.Age;
swarmchart(FC,Age,'filled','CData',ManuscriptColors.CutotypeColors(FC,:))
xlim([0 4])
hold on
arrayfun(@(x) plot([x-.4 x+.4],[median(Age(FC==x)) median(Age(FC==x))],'Color','black','LineWidth',1), 1:3)
pbaspect([1 1 1])
xticks(1:3)
xticklabels(["FC1-Children" "FC2-children" "Parents"])
ylabel('Age at sampling (years)')
[p12] = ranksum(Age(FC==1),Age(FC==2));
[p13] = ranksum(Age(FC==1),Age(FC==3));
[p23] = ranksum(Age(FC==2),Age(FC==3));

text(1.5,65,string(p12))
text(2.5,65,string(p23))
text(2,70,string(p13))
