function M = BrackenDataUMAP(M,ManuscriptColors)
% This module takes M, the subject-level metagenomics sample table, and 
% clusters samples with UMAP at the species-level in addition to various 
% supplementary figures 

% for UMAP
path('./umap',path);
path('./epp',path);


%% Get the ages, and colors based off of the ages



Age=M.Age;
GroupColors=getGroupProgressiveColors(Age,ManuscriptColors);



%% Get species-level bracken abundances, and remove sparse taxa and human-associated taxa



Names=M.NamesSpecies(1,:);
A = M.BrackenSpecies;    

% remove human associated taxa

HumanTaxa = contains(Names,"Human")|contains(Names,"Toxoplasma");
HighestAssigned = max(A,[],1);
Remove = HumanTaxa|HighestAssigned<.005;
A = A(:,~Remove);
Names=Names(~Remove);
A=A./sum(A,2);



%% run UMAP



figure
[UMAPCoordinates,UMAPParams,UMAPClusters,ClusterExtras]=run_umap(A,'min_dist',.1,'init','spectral','n_epochs',10000,'n_components',2,'metric','chebychev','cluster_detail','adaptive','method','MEX');



%% add UMAP data to M



M.BrackenAbundanceUMAPClustering=A;
M.BrackenNamesUMAPClustering=repmat(Names,size(M,1),1);



%% Save UMAP data



save('data/UMAPrunInfo.mat','UMAPCoordinates','UMAPParams','UMAPClusters','ClusterExtras')



%% Add clusters to M



M.UMAPClustersSubjectTime=UMAPClusters';



%% make cutotypes for each subject-- subject is FC2 if at least one timepoint is FC2



UMAPClustersSubject=arrayfun(@(x) {M.UMAPClustersSubjectTime(x==M.SID)}, M.SID);
UMAPClustersSubject=cellfun(@max, UMAPClustersSubject);
M.UMAPClustersSubject=UMAPClustersSubject;



%% plot UMAP colored by age, shapes are sex

Male = M.Sex=="M";

f=figure;
hold on
scatter(UMAPCoordinates(Male,1),UMAPCoordinates(Male,2),'filled','Cdata',GroupColors(Male,:),'Marker','^')
scatter(UMAPCoordinates(~Male,1),UMAPCoordinates(~Male,2),'filled','Cdata',GroupColors(~Male,:),'Marker','o')
pbaspect([1 1 1])
xlabel('UMAP 1')
ylabel('UMAP 2')

saveas(f,'ManuscriptFigures/SpeciesLevelSubjectTimeUMAP_colored_by_age.svg')



%% get species and genus-level data of the three main taxa



AbundanceCacnes = A(:,Names=="Cutibacterium acnes");
AbundanceSepi = A(:,Names=="Staphylococcus epidermidis");
AbundanceSmitis= A(:,Names=="Streptococcus mitis");
AbundanceStrep = sum(A(:,startsWith(Names,"Streptococcus")),2);
AbundanceStaph = sum(A(:,startsWith(Names,"Staphylococcus")),2);
AbundanceCuti = sum(A(:,startsWith(Names,"Cutibacterium")),2);

RelaventTaxaNames = ["Cutibacterium acnes" "Staphylococcus epidermidis" "Streptococcus mitis" "Cutibacterium" "Staphylococcus" "Streptococcus"];
RelaventTaxaAbundances = [AbundanceCacnes AbundanceSepi AbundanceSmitis AbundanceCuti AbundanceStaph AbundanceStrep];

%%


MeanAbundanceCacnes=arrayfun(@(x) mean(AbundanceCuti(UMAPClusters==x)), [1 2])
%%
MeanAbundanceSepi=arrayfun(@(x) mean(AbundanceSepi(UMAPClusters==x)), [1 2])

%% UMAP colored by abundance of different taxa



f=figure;

Male = M.Sex=="M"
for i =1:numel(RelaventTaxaNames)
    Ms = find(Male);
    Fs = find(~Male);
    [~,Midx] = sort(RelaventTaxaAbundances(Ms,i),'ascend');
    [~,Fidx] = sort(RelaventTaxaAbundances(Fs,i),'ascend');
    Ms = Ms(Midx);
    Fs=Fs(Fidx);
    subplot(2,3,i)
    clrs=getLinearColors(RelaventTaxaAbundances(:,i),'Greens',20);
    s=scatter(UMAPCoordinates(Ms,1),UMAPCoordinates(Ms,2),18,'filled','CData',clrs(Ms,:),'Marker','^');
    hold on
    [~,pv,~]=kstest2(RelaventTaxaAbundances(UMAPClusters==1,i),RelaventTaxaAbundances(UMAPClusters==2,i))
    s=scatter(UMAPCoordinates(Fs,1),UMAPCoordinates(Fs,2),18,'filled','CData',clrs(Fs,:),'Marker','o');
    title([char(RelaventTaxaNames(i))])
    xlabel('UMAP1');ylabel('UMAP2');
    text(32,-8,string(pv))
end

%

saveas(f,'ManuscriptFigures/SpeciesLevelSubjectTimeUMAP_colored_by_taxa.svg')



%% histograms of taxa abundance between clusters



f=figure;
c=1
bw=.05
for i =1:numel(RelaventTaxaNames)
    
    
    Abundances = RelaventTaxaAbundances(:,i);
    subplot(2,3,i)
    hold on
    histogram(Abundances(M.UMAPClustersSubjectTime==1),'FaceColor',ManuscriptColors.CutotypeColors(2,:),'BinWidth',bw,'Normalization','probability')
    histogram(Abundances(M.UMAPClustersSubjectTime==2),'FaceColor',ManuscriptColors.CutotypeColors(3,:),'BinWidth',bw,'Normalization','probability')
    xlim([0 1])
    ylim([0 1])
    yticks([0:.5:1])
    title([char(RelaventTaxaNames(i))])
    pbaspect([1 1 1])
    if i==1|i==4
        ylabel('proportion samples')
        l=legend
        l.String={'FC1' 'FC2'}
    end
end

saveas(f,'ManuscriptFigures/taxa_abundance_by_cutotype.svg')



%% age histograms 

uSID = unique(M.SID);
Cutotype = arrayfun(@(x) unique(M.subject_cutotype(M.SID==x)), uSID);
AverageAge= arrayfun(@(x) mean(M.Age(M.SID==x)), uSID);

f=figure
subplot(2,1,1)
histogram(AverageAge(Cutotype==1),'BinWidth',2,'FaceColor',[.85 .85 .85],'Normalization','probability')
xlim([0 65])
ylim([0 .4])
subplot(2,1,2)
histogram(AverageAge(Cutotype==2),'BinWidth',2,'FaceColor',[.85 .85 .85],'Normalization','probability')
xlim([0 65])
ylim([0 .4])
%
saveas(f,'ManuscriptFigures/AgeHistograms.svg')



%% Ages of child subjects only 



f=figure;
Child = AverageAge<18;
clrs = ManuscriptColors.CutotypeColors;
clrs = clrs(Cutotype,:);
swarmchart(Cutotype(Child),AverageAge(Child),'filled','CData',clrs(Child,:),'XJitter','density',XJitterWidth=.5)
xlim([.5 2.5])
ylim([0 15])
yticks([0:5:15])
xticks([1:2])
xticklabels(["FC1 Children" "FC2 Children"])
pv = ranksum(AverageAge(Child&Cutotype==1),AverageAge(Child&Cutotype==2));
text(1.5, 15,['p=' char(string(round(pv,2,'significant')))])
hold on
arrayfun(@(x) plot([x-.25 x+.25], [median(AverageAge(Child&Cutotype==x)) median(AverageAge(Child&Cutotype==x))],'Color','black','LineWidth',1)   , 1:2)
%
saveas(f,'ManuscriptFigures/AgeSwarmchartsChildrenOnly.svg')


%% plot average abundance of relavent taxa against relavent ages



AverageAbundancesRelaventTaxa = arrayfun(@(x) {mean(RelaventTaxaAbundances(M.SID==x,:),1)}, uSID);
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

%%



saveas(h,'ManuscriptFigures/AbundanceVersusAge.svg')



%% Bi-plot of Species-level abundance data

%Species-level PCA 
[coeff,score,~,~,explained,~] = pca(A,"NumComponents",2);

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
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Names(Biggest3),'Color','black','LineWidth',1);
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

saveas(f,'ManuscriptFigures/Figure_S5_A_BiPlotSpeciesLevel','svg')

%% As above, but genus level


Names=M.NamesGenus(1,:);
A = M.BrackenGenus;    

% remove human associated taxa

HumanTaxa = contains(Names,"Human")|contains(Names,"Toxoplasma");
HighestAssigned = max(A,[],1);
Remove = HumanTaxa|HighestAssigned<.005;
A = A(:,~Remove);
Names=Names(~Remove);
A=A./sum(A,2);




%Species-level PCA 
[coeff,score,~,~,explained,~] = pca(A,"NumComponents",2);

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
biplot(coeff(Biggest3,[1 2]),'MarkerSize',5,'VarLabels',Names(Biggest3),'Color','black','LineWidth',1);
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



saveas(f,'ManuscriptFigures/Figure_S5_B_BiPlotGenusLevel','svg')

%%


