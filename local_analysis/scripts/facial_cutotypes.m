function MsubjectTime = facial_cutotypes(MsubjectTime,ManuscriptColors)
% This module takes MsubjectTime, the subject-level metagenomics sample table, and
% clusters samples at the species-level


%% Get the ages, and colors based off of the ages

GroupColors=getGroupProgressiveColors(MsubjectTime.Age,ManuscriptColors);
Male = MsubjectTime.Sex=="M"; % boolean for plotting

%% Get species-level bracken abundances, and remove sparse taxa and human-associated taxa

SpeciesNames=MsubjectTime.NamesSpecies(1,:);
SpeciesAbundance = MsubjectTime.BrackenSpecies;

% remove human associated taxa
% Homo sapiens DNA is not here, but there is some other human associated
% viral things here.
HumanTaxa = contains(SpeciesNames,"Human")|contains(SpeciesNames,"Toxoplasma");
HighestAssigned = max(SpeciesAbundance,[],1);
Remove = HumanTaxa|HighestAssigned<.005;
SpeciesAbundance = SpeciesAbundance(:,~Remove);
SpeciesNames=SpeciesNames(~Remove);
SpeciesAbundance=SpeciesAbundance./sum(SpeciesAbundance,2);

%% Get Euclidian pairwise distances amongst subjects

D=(pdist(SpeciesAbundance,'euclidean'));

%% Measure optimal number of clusters using different criteria

XX=50;
e_sil=evalclusters(SpeciesAbundance,'linkage','silhouette','KList',1:XX,'Distance','Euclidean');
e_DaviesBouldin=evalclusters(SpeciesAbundance,'linkage','DaviesBouldin','KList',1:XX);
e_CalinskiHarabasz=evalclusters(SpeciesAbundance,'linkage','CalinskiHarabasz','KList',1:XX);

f=figure;
subplot(3,1,1)
scatter(1:XX,e_sil.CriterionValues,'filled','MarkerFaceColor','black')
best = max(e_sil.CriterionValues);
hold on
scatter(find(e_sil.CriterionValues==best),e_sil.CriterionValues(e_sil.CriterionValues==best),'filled','MarkerFaceColor','r')
xlabel('K clusters')
xlim([0 50])
xticks(0:10:50)
ylabel('Silhoette metric')
pbaspect([1 1 1])
subplot(3,1,2)
scatter(1:XX,e_DaviesBouldin.CriterionValues,'filled','MarkerFaceColor','black')
best = min(e_DaviesBouldin.CriterionValues);
hold on
scatter(find(e_DaviesBouldin.CriterionValues==best),e_DaviesBouldin.CriterionValues(e_DaviesBouldin.CriterionValues==best),'filled','MarkerFaceColor','r')
xlabel('K clusters')
ylabel('Davies Bouldin metric')
xlim([0 50])
xticks(0:10:50)
pbaspect([1 1 1])
subplot(3,1,3)
scatter(1:XX,e_CalinskiHarabasz.CriterionValues,'filled','MarkerFaceColor','black')
best = max(e_CalinskiHarabasz.CriterionValues);
hold on
scatter(find(e_CalinskiHarabasz.CriterionValues==best),e_CalinskiHarabasz.CriterionValues(e_CalinskiHarabasz.CriterionValues==best),'filled','MarkerFaceColor','r')
xlabel('K clusters')
ylabel('Calinski Harabasz metric')
pbaspect([1 1 1])
xlim([0 50])
xticks(0:10:50)

yticks([0:50:250])

%% Show that both types of linkage make same clusters

Lcomplete=linkage(D,'complete');
clusters_completelinkage = cluster(Lcomplete,'MaxClust',2);
Lward=linkage(D,'ward');
clusters_wardlinkage = cluster(Lward,'MaxClust',2);

% expression is true where clusters are the same or inverse, i.e. the same
% but with different numbering (numbers are arbitrary)

clusters_are_the_same = all(clusters_wardlinkage==clusters_completelinkage)|(all(clusters_wardlinkage~=clusters_completelinkage)&all(clusters_completelinkage==1|clusters_completelinkage==2)&all(clusters_completelinkage==1|clusters_completelinkage==2));
if clusters_are_the_same
    fprintf('Both complete and ward linkage result in the same cluster membership \n')
end

f=figure;
subplot(2,1,1)
d=dendrogram(Lcomplete,Inf,'Labels',string(clusters_completelinkage));
xlabel('complete linkage tree of each of 101 clustered points')
pbaspect([1 1 1])
subplot(2,1,2)
d=dendrogram(Lward,Inf,'Labels',string(clusters_wardlinkage));
xlabel('ward linkage tree of each of 101 clustered points')
pbaspect([1 1 1])

%% Get facial cutotype clusters

clusters_subjecttime_bracken_species = clusters_wardlinkage;

%% add UMAP data to M


MsubjectTime.BrackenAbundanceClustering=SpeciesAbundance;
MsubjectTime.BrackenNamesClustering=repmat(SpeciesNames,size(MsubjectTime,1),1);


%% Add clusters to M


MsubjectTime.ClustersSubjectTime=clusters_subjecttime_bracken_species;


%% make cutotypes for each subject-- subject is FC2 if at least one timepoint is FC2


ClustersSubject=arrayfun(@(x) {MsubjectTime.ClustersSubjectTime(x==MsubjectTime.SID)}, MsubjectTime.SID);
ClustersSubject=cellfun(@max, ClustersSubject);
MsubjectTime.ClustersSubject=ClustersSubject;
MsubjectTime.clusters_subjecttime_bracken_species=MsubjectTime.ClustersSubjectTime;


%% get species and genus-level data of the three main taxa


AbundanceCacnes = SpeciesAbundance(:,SpeciesNames=="Cutibacterium acnes");
AbundanceSepi = SpeciesAbundance(:,SpeciesNames=="Staphylococcus epidermidis");
AbundanceSmitis= SpeciesAbundance(:,SpeciesNames=="Streptococcus mitis");
AbundanceStrep = sum(SpeciesAbundance(:,startsWith(SpeciesNames,"Streptococcus")),2);
AbundanceStaph = sum(SpeciesAbundance(:,startsWith(SpeciesNames,"Staphylococcus")),2);
AbundanceCuti = sum(SpeciesAbundance(:,startsWith(SpeciesNames,"Cutibacterium")),2);

RelaventTaxaNames = ["Cutibacterium acnes" "Staphylococcus epidermidis" "Streptococcus mitis" "Cutibacterium" "Staphylococcus" "Streptococcus"];
RelaventTaxaAbundances = [AbundanceCacnes AbundanceSepi AbundanceSmitis AbundanceCuti AbundanceStaph AbundanceStrep];

%% Ages of child subjects only

[uSID,iSID,~]= unique(MsubjectTime.SID);
Cutotype = arrayfun(@(x) unique(MsubjectTime.ClustersSubject(MsubjectTime.SID==x)), uSID);
AverageAge= arrayfun(@(x) mean(MsubjectTime.Age(MsubjectTime.SID==x)), uSID);
male_subject = Male(iSID);

f=figure;
Child = AverageAge<18;
clrs = ManuscriptColors.CutotypeColors;
clrs = clrs(Cutotype,:);

s=swarmchart(Cutotype(Child&male_subject),AverageAge(Child&male_subject),'filled','CData',clrs(Child&male_subject,:),'Marker','^','XJitter','density',XJitterWidth=.5)
hold on
s=swarmchart(Cutotype(Child&~male_subject),AverageAge(Child&~male_subject),'filled','CData',clrs(Child&~male_subject,:),'Marker','o','XJitter','density',XJitterWidth=.5)

xlim([.5 2.5])
ylim([0 15])
yticks([0:5:15])
xticks([1:2])
xticklabels(["FC1 Children" "FC2 Children"])
pv = ranksum(AverageAge(Child&Cutotype==1),AverageAge(Child&Cutotype==2));

text(1.5, 15,['p=' char(string(round(pv,2,'significant')))])
hold on
arrayfun(@(x) plot([x-.25 x+.25], [median(AverageAge(Child&Cutotype==x)) median(AverageAge(Child&Cutotype==x))],'Color','black','LineWidth',1)   , 1:2)

%% plot average abundance of relavent taxa against relavent ages


AverageAbundancesRelaventTaxa = arrayfun(@(x) {mean(RelaventTaxaAbundances(MsubjectTime.SID==x,:),1)}, uSID);
AverageAbundancesRelaventTaxa=vertcat(AverageAbundancesRelaventTaxa{:});

AbundanceCuti=AverageAbundancesRelaventTaxa(:,1);
AbundanceSepi=AverageAbundancesRelaventTaxa(:,2);
AbundanceStrep=AverageAbundancesRelaventTaxa(:,6);

Cutotype(AverageAge>18)=3;

clrs3 = ManuscriptColors.CutotypeColors;
clrs3=clrs3(Cutotype,:)

h=figure
X = AverageAge;
bool = X<20;
%
subplot(3,2,1)
hold on
scatter(X(bool),AbundanceCuti(bool),'filled','CData',clrs3(bool,:))
[f,gof] = fit(X(bool),AbundanceCuti(bool),'exp1');
plot(f)
text(5,1,['adjR2=' string(gof.adjrsquare)])
[~,pv]=corr(X(bool),AbundanceCuti(bool))
text(5,.8,string(pv));
xlim([5 15])
ylim([0 1])
pbaspect([1 1 1])
%
subplot(3,2,3)
scatter(X(bool),AbundanceSepi(bool),'filled','CData',clrs3(bool,:))
[f,gof] = fit(X(bool),AbundanceSepi(bool),'exp1')
hold on
plot(f)
xlim([5 15])
ylim([0 .5])
pbaspect([1 1 1])
text(5,.5,['adjR2=' string(gof.adjrsquare)])
[~,pv]=corr(X(bool),AbundanceSepi(bool))
text(5,.4,string(pv));
yticks(0:.25:.5)
subplot(3,2,5)
scatter(X(bool),AbundanceStrep(bool),'filled','CData',clrs3(bool,:))
xlim([5 15])
ylim([0 .5])
pbaspect([1 1 1])
[f,gof] = fit(X(bool),AbundanceStrep(bool),'poly1')
hold on
plot(f)
text(5,.5,['adjR2=' string(gof.adjrsquare)])
yticks(0:.25:.5)
[~,pv]=corr(X(bool),AbundanceStrep(bool))
text(5,.4,string(pv));%
%
X = AverageAge;
bool = X>20;
subplot(3,2,2)
scatter(X(bool),AbundanceCuti(bool),'filled','CData',clrs3(bool,:))
xlim([32 62])
ylim([0 1])
pbaspect([1 1 1])
[f,gof] = fit(X(bool),AbundanceCuti(bool),'exp1')
hold on
plot(f)
text(5,1,['adjR2=' string(gof.adjrsquare)])
[~,pv]=corr(X(bool),AbundanceCuti(bool))
text(5,.8,string(pv));
%
subplot(3,2,4)
scatter(X(bool),AbundanceSepi(bool),'filled','CData',clrs3(bool,:))
xlim([32 62])
ylim([0 .5])
hold on
[f,gof] = fit(X(bool),AbundanceSepi(bool),'exp1')
plot(f)
pbaspect([1 1 1])
text(5,.5,['adjR2=' string(gof.adjrsquare)])
yticks(0:.25:.5)
[~,pv]=corr(X(bool),AbundanceStaph(bool))
text(5,.4,string(pv));
subplot(3,2,6)
scatter(X(bool),AbundanceStrep(bool),'filled','CData',clrs3(bool,:))
xlim([32 62])
ylim([0 .5])
yticks(0:.25:.5)
pbaspect([1 1 1])
hold on
[f,gof] = fit(X(bool),AbundanceStrep(bool),'exp1')
plot(f)
text(5,.5,['adjR2=' string(gof.adjrsquare)])
[~,pv]=corr(X(bool),AbundanceStrep(bool))
text(5,.4,string(pv));

%% Species-level PCA

[coeff,score,~,~,explained,~] = pca(SpeciesAbundance,"NumComponents",2);

%% plot just PCA for main text

f=figure;
% makes the scatter plot from score
clrs = GroupColors;
sc1=scatter(score(Male,1),score(Male,2),18,'filled','CData',clrs(Male,:),'LineWidth',.25,'Marker','^');
hold on
sc2=scatter(score(~Male,1),score(~Male,2),18,'filled','CData',clrs(~Male,:),'LineWidth',.25,'Marker','o');
% Choose which taxa you want to label on plot (too many to have them all)
[~,ia]=sort(sqrt((coeff(:,1).^2+coeff(:,2).^2)),'descend');
Biggest3 = ia(1:3);
% First, plot the ones which are not the biggest 3
idx=1:size(coeff,1);
% biplot(coeff(~ismember(idx,Biggest3),[1 2]),'MarkerSize',5,'VarLabels',"",'Color','black','LineWidth',1);
% % make the second part, which is the scores for each sample. color by age group
% % Then, plot the 3 most important taxa with the label
% % biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Taxa(Biggest3),'Color','black','LineWidth',1);
% biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Names(Biggest3),'Color','black','LineWidth',1);
% add labels and set axis limits
xlim([-.5 .5])
ylim([-.5 .5])
h=gca;
h.LineWidth=1;
xticks([-1:.5:1])
yticks([-1:.5:1])
pbaspect([1 1 1])
xlabel({['Component 1(' char(string(round(explained(1),2))) ')%'] ; 'Colored by age'})
ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])

%% Plot species-level PCA

f=figure;

subplot(3,1,1)
% makes the scatter plot from score
clrs = GroupColors;
sc1=scatter(score(Male,1),score(Male,2),18,'filled','CData',clrs(Male,:),'LineWidth',.25,'Marker','^');
hold on
sc2=scatter(score(~Male,1),score(~Male,2),18,'filled','CData',clrs(~Male,:),'LineWidth',.25,'Marker','o');
% Choose which taxa you want to label on plot (too many to have them all)
[~,ia]=sort(sqrt((coeff(:,1).^2+coeff(:,2).^2)),'descend');
Biggest3 = ia(1:3);
% First, plot the ones which are not the biggest 3
idx=1:size(coeff,1);
biplot(coeff(~ismember(idx,Biggest3),[1 2]),'MarkerSize',5,'VarLabels',"",'Color','black','LineWidth',1);
% make the second part, which is the scores for each sample. color by age group
% Then, plot the 3 most important taxa with the label
% biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Taxa(Biggest3),'Color','black','LineWidth',1);
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',SpeciesNames(Biggest3),'Color','black','LineWidth',1);
% add labels and set axis limits
xlim([-1 1])
ylim([-1 1])
h=gca;
h.LineWidth=1;
xticks([-1:.5:1])
yticks([-1:.5:1])
pbaspect([1 1 1])
xlabel({['Component 1(' char(string(round(explained(1),2))) ')%'] ; 'Colored by age'})
ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])
%
subplot(3,1,2)
clrs = zeros(101,3);
clrs(MsubjectTime.ClustersSubjectTime==1,:)=.9;
sc1=scatter(score(Male,1),score(Male,2),18,'filled','CData',clrs(Male,:),'LineWidth',.25,'Marker','^');
hold on
sc2=scatter(score(~Male,1),score(~Male,2),18,'filled','CData',clrs(~Male,:),'LineWidth',.25,'Marker','o');
hold on;
% Choose which taxa you want to label on plot (too many to have them all)
[~,ia]=sort(sqrt((coeff(:,1).^2+coeff(:,2).^2)),'descend');
Biggest3 = ia(1:3);
% First, plot the ones which are not the biggest 3
idx=1:size(coeff,1);
biplot(coeff(~ismember(idx,Biggest3),[1 2]),'MarkerSize',5,'VarLabels',"",'Color','black','LineWidth',1);
% make the second part, which is the scores for each sample. color by age group
% Then, plot the 3 most important taxa with the label
% biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Taxa(Biggest3),'Color','black','LineWidth',1);
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',SpeciesNames(Biggest3),'Color','black','LineWidth',1);
% add labels and set axis limits
xlim([-1 1])
ylim([-1 1])
h=gca;
h.LineWidth=1;
xticks([-1:.5:1])
yticks([-1:.5:1])
pbaspect([1 1 1])
xlabel({['Component 1(' char(string(round(explained(1),2))) ')%'] ; 'Colored by FC'})
ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])

%
subplot(3,1,3)
clrs = zeros(101,3);
clrs(Male,:)=.9;
sc1=scatter(score(Male,1),score(Male,2),18,'filled','CData',clrs(Male,:),'LineWidth',.25,'Marker','^');
hold on
sc2=scatter(score(~Male,1),score(~Male,2),18,'filled','CData',clrs(~Male,:),'LineWidth',.25,'Marker','o');
hold on;
% Choose which taxa you want to label on plot (too many to have them all)
[~,ia]=sort(sqrt((coeff(:,1).^2+coeff(:,2).^2)),'descend');
Biggest3 = ia(1:3);
% First, plot the ones which are not the biggest 3
idx=1:size(coeff,1);
biplot(coeff(~ismember(idx,Biggest3),[1 2]),'MarkerSize',5,'VarLabels',"",'Color','black','LineWidth',1);
% make the second part, which is the scores for each sample. color by age group
% Then, plot the 3 most important taxa with the label
% biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Taxa(Biggest3),'Color','black','LineWidth',1);
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',SpeciesNames(Biggest3),'Color','black','LineWidth',1);
% add labels and set axis limits
xlim([-1 1])
ylim([-1 1])
h=gca;
h.LineWidth=1;
xticks([-1:.5:1])
yticks([-1:.5:1])
pbaspect([1 1 1])
xlabel({['Component 1(' char(string(round(explained(1),2))) ')%'] ; 'Colored by Sex'})
ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])

%% species PCA colored by different species

f=figure;

clusters=MsubjectTime.clusters_subjecttime_bracken_species;
for i =1:numel(RelaventTaxaNames)
    Ms = find(Male);
    Fs = find(~Male);
    [~,Midx] = sort(RelaventTaxaAbundances(Ms,i),'ascend');
    [~,Fidx] = sort(RelaventTaxaAbundances(Fs,i),'ascend');
    Ms = Ms(Midx);
    Fs=Fs(Fidx);
    subplot(2,3,i)
    clrs=getLinearColors(RelaventTaxaAbundances(:,i),'Greens',20);
    s=scatter(score(Ms,1),score(Ms,2),18,'filled','CData',clrs(Ms,:),'Marker','^');
    hold on
    pv=ranksum(RelaventTaxaAbundances(clusters==1,i),RelaventTaxaAbundances(clusters==2,i));
    ['Median FC1 abundance of ' char(RelaventTaxaNames(i)) ' = ' char(string(mean(RelaventTaxaAbundances(clusters==1,i))))]
    ['Median FC2 abundance of ' char(RelaventTaxaNames(i)) ' = ' char(string(mean(RelaventTaxaAbundances(clusters==2,i))))]

    s=scatter(score(Fs,1),score(Fs,2),18,'filled','CData',clrs(Fs,:),'Marker','o');
    title([char(RelaventTaxaNames(i))])
    xlim([-1 1])
    ylim([-1 1])
    h=gca;
    h.LineWidth=1;
    xticks([-1:.5:1])
    yticks([-1:.5:1])
    pbaspect([1 1 1])
    xlabel(['Component 1(' char(string(round(explained(1),2))) ')%'])
    ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])
    text(0,.9,['p=' string(pv)])
end

%% As above, but genus level


GenusNames=MsubjectTime.NamesGenus(1,:);
GenusAbundance = MsubjectTime.BrackenGenus;

% remove human associated taxa

HumanTaxa = contains(GenusNames,"Human")|contains(GenusNames,"Toxoplasma");
HighestAssigned = max(GenusAbundance,[],1);
Remove = HumanTaxa|HighestAssigned<.005;
GenusAbundance = GenusAbundance(:,~Remove);
GenusNames=GenusNames(~Remove);
GenusAbundance=GenusAbundance./sum(GenusAbundance,2);
%Species-level PCA
[coeff,score,~,~,explained,~] = pca(GenusAbundance,"NumComponents",2);
%
% Initialize a figure
f=figure;

% makes the scatter plot from score
sc=scatter(score(:,1),score(:,2),'filled','CData',GroupColors,'LineWidth',.25);
sc.SizeData=18;
hold on;
% Choose which taxa you want to label on plot (too many to have them all)
[~,ia]=sort(sqrt((coeff(:,1).^2+coeff(:,2).^2)),'descend');
Biggest3 = ia(1:3);
% First, plot the ones which are not the biggest 3
idx=1:size(coeff,1);
biplot(coeff(~ismember(idx,Biggest3),[1 2]),'MarkerSize',5,'VarLabels',"",'Color','black','LineWidth',1);
% make the second part, which is the scores for each sample. color by age group
% Then, plot the 3 most important taxa with the label
% biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Taxa(Biggest3),'Color','black','LineWidth',1);
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',GenusNames(Biggest3),'Color','black','LineWidth',1);
% add labels and set axis limits
xlim([-1 1])
ylim([-1 1])
h=gca;
h.LineWidth=1;
xticks([-1:.5:1])
yticks([-1:.5:1])
pbaspect([1 1 1])
xlabel(['Component 1(' char(string(round(explained(1),2))) ')%'])
ylabel(['Component 2(' char(string(round(explained(2),2))) ')%'])

end
