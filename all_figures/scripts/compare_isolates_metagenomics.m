function compare_isolates_metagenomics(IT,M,cutoff_isolates,cutoff_assigned,cutoff_real)
%% Get abundances to use in comparison for both species



MetagenomicsCacnesAbundance = M.CombinedCacnesLineages;
MetagenomicsSepiAbundance=M.CombinedSepiLineages;



%% get cluster numbers



MetagenomicsClustersCacnes=M.CacnesLineageNumbers;
MetagenomicsClustersSepi=M.SepiLineageNumbers;



%% comparing subject-time (samplings), so get this for each row in M



MetagenomicsSubjectTime=arrayfun(@(x) strjoin([M.SID(x) string(M.TP(x))],'') , 1:size(M,1))';



%% run function to get numbers



CacnesIsolateNumber=mg2isolates(MetagenomicsClustersCacnes,MetagenomicsSubjectTime,IT(IT.SpeciesString=="cacnes",:));
SepiIsolateNumber=mg2isolates(MetagenomicsClustersSepi,MetagenomicsSubjectTime,IT(IT.SpeciesString=="sepi",:));



%% cutoffs



%%
[p1,p2, p3, p4]= plot_comparisons(CacnesIsolateNumber,MetagenomicsCacnesAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.85 .85 .85]);
%


saveas(p1,'ManuscriptFigures/CacnesNlineagesMG_Iso.svg')
saveas(p2,'ManuscriptFigures/supplemental_CacnesAbundanceNisoVsAbundanceMG.svg')
saveas(p3,'ManuscriptFigures/supplemental_CacnesAbundanceUnfoundLineages.svg')
% saveas(p4,'ManuscriptFigures/supplemental_CacnesAbundanceExcludedLineages.svg')



%%


[p1,p2, p3, p4]= plot_comparisons(SepiIsolateNumber,MetagenomicsSepiAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.25 .25 .25]);

saveas(p1,'ManuscriptFigures/SepiNlineagesMG_Iso.svg')
saveas(p2,'ManuscriptFigures/supplemental_SepiAbundanceNisoVsAbundanceMG.svg')
saveas(p3,'ManuscriptFigures/supplemental_SepiAbundanceUnfoundLineages.svg')
% saveas(p4,'ManuscriptFigures/supplemental_SepiAbundanceExcludedLineages.svg')

%%

p1 = plot_both(CacnesIsolateNumber,MetagenomicsCacnesAbundance,SepiIsolateNumber,MetagenomicsSepiAbundance,cutoff_isolates,cutoff_assigned,cutoff_real)

%%
p2 = plot_both_histograms(CacnesIsolateNumber,MetagenomicsCacnesAbundance,SepiIsolateNumber,MetagenomicsSepiAbundance,cutoff_isolates,cutoff_assigned,cutoff_real)

%%
saveas(p1,'ManuscriptFigures/NlineagesMetagenomicsIsolates.svg')
saveas(p2,'ManuscriptFigures/NLineagesHistogram.svg')

end

%% Get Isolate numbers array corrosponding to metagenomics abundances

function IsolateNumber=mg2isolates(Clusters,SubjectTime,IT)

% size of array to be made (isolate output same size as one from MG)
[I,J]=size(Clusters);
N=I*J;

SubjectTime=repmat(SubjectTime,1,J);


IsolateNumber=arrayfun(@(x) sum(IT.SubjectTime==SubjectTime(x)&IT.ClusterString==Clusters(x)), 1:N);

IsolateNumber=reshape(IsolateNumber,I,J);

end

%% make the scatter plot


function [p1,p2, p3,p4]= plot_comparisons(Isolates,Metagenomics,cutoff_isolates,cutoff_assigned,cutoff_real,clr)

%% remove rows with insufficient data

AnyIsolatesAtAll=sum(Isolates,1)>0;
MG_excluded_lineages=Metagenomics(:,~AnyIsolatesAtAll);

Isolates=Isolates(:,AnyIsolatesAtAll);
Metagenomics=Metagenomics(:,AnyIsolatesAtAll);

EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);



%%



Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);

% scatter plot
N = numel(Niso);
jitter = rand(N,2)/3;

%% Number of lineages in isolates vs metagenomics
p1=figure
Xm=10;
scatter(Niso+jitter(:,1),NMetagenomics+jitter(:,2),'filled','MarkerFaceColor',clr);
hold on
plot([0 Xm],[0 Xm],'Color','black')
xlim([0 Xm])
ylim([0 Xm])
pbaspect([1 1 1])
[rho,p]=corr(Niso,NMetagenomics,'type','Pearson')
text(0,Xm,['p=' char(string(p))])
text(0,Xm-2,['R2=' char(string(rho^2))])

xlabel('Number of lineages (Isolates)')
ylabel('Number of lineages (Metagenomics)')

%% Abundance of lineages found at both tps

FoundBoth=(Isolates>0&Metagenomics>0);
MetagenomicsBoth=Metagenomics(FoundBoth);
IsolatesBoth=Isolates(FoundBoth);
p2=figure;

scatter(IsolatesBoth,MetagenomicsBoth,'filled')
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
xlabel("Relative Abundance (isolates)")
ylabel("Relative Abundance (metagenomics)")
title('lineages found in both datasets')
hold on
plot([0 1],[0 1],'Color','black')

[rho,p]=corr(IsolatesBoth,MetagenomicsBoth,'type','Pearson')
text(0,1,['p=' char(string(p))])
text(0,.8,['R2=' char(string(rho^2))])

%% abundance of lineages not found in other datasets

p3=figure

MetagenomicsMissingInIso=Metagenomics(Metagenomics>0&Isolates==0);
IsoMissingInMetagenomics=Isolates(Isolates>0&Metagenomics==0);

subplot(2,1,1)
histogram(MetagenomicsMissingInIso,'BinWidth',.01)
xlim([0 1])
% ylim([0 1])
pbaspect([1 1 1])
xlabel('Relative abundance (Metagenomics) of lineages not found in isolates and found in metagenomics')
ylabel('number of instances')
xticks([0:.2:1])
subplot(2,1,2)

histogram(IsoMissingInMetagenomics,'BinWidth',.1)
xlim([0 1])
% ylim([0 1])
xticks([0:.2:1])
pbaspect([1 1 1])
xlabel('Relative abundance (isolates) of lineages not found in metagenomcs and found in isolates')
ylabel('number of instances')

%% lneages which were exlcued

p4=figure;
histogram(MG_excluded_lineages(MG_excluded_lineages>0),'BinWidth',.01)
xlabel('Nonzero abundances of excluded lineages')
ylabel('number of instances')
pbaspect([1 1 1])
xlim([0 1])
% ylim([0 1])
end

%%


function p1=plot_both(Isolates1,Metagenomics1,Isolates2,Metagenomics2,cutoff_isolates,cutoff_assigned,cutoff_real)
p1=figure

Isolates=Isolates1;
Metagenomics=Metagenomics1;

AnyIsolatesAtAll=sum(Isolates,1)>0;

Isolates=Isolates(:,AnyIsolatesAtAll);
Metagenomics=Metagenomics(:,AnyIsolatesAtAll);

EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);
%%
Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);
found_either = min([Niso NMetagenomics],[],2)
median_min_found=median(found_either)
% scatter plot
N = numel(Niso);
jitter = rand(N,2)/3;
% Number of lineages in isolates vs metagenomics
Xm=10;
scatter(Niso+jitter(:,1),NMetagenomics+jitter(:,2),'filled','MarkerFaceColor',[.85 .85 .85]);
[rho,p]=corr(Niso,NMetagenomics,'type','Pearson')
text(0,Xm,['p_cacnes_gray=' char(string(p))])
text(0,Xm-1,['R2_cacnes_gray=' char(string(rho^2))])
hold on

%%
Isolates=Isolates2;
Metagenomics=Metagenomics2;

AnyIsolatesAtAll=sum(Isolates,1)>0;
MG_excluded_lineages=Metagenomics(:,~AnyIsolatesAtAll);

Isolates=Isolates(:,AnyIsolatesAtAll);
Metagenomics=Metagenomics(:,AnyIsolatesAtAll);

EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);
%%
Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);

% scatter plot
N = numel(Niso);
jitter = rand(N,2)/3;
% Number of lineages in isolates vs metagenomics
Xm=10;
scatter(Niso+jitter(:,1),NMetagenomics+jitter(:,2),'filled','MarkerFaceColor',[0 0 0]);
hold on
%%
plot([0 Xm],[0 Xm],'Color','black')
xlim([0 Xm])
ylim([0 Xm])
pbaspect([1 1 1])
[rho,p]=corr(Niso,NMetagenomics,'type','Pearson')
text(0,Xm-2,['p_sepi_black=' char(string(p))])
text(0,Xm-3,['R2_sepi_black=' char(string(rho^2))])

found_either = min([Niso NMetagenomics],[],2)
median_min_found=median(found_either)
xlabel('Number of lineages (Isolates)')
ylabel('Number of lineages (Metagenomics)')


end

%%



function p1=plot_both_histograms(Isolates1,Metagenomics1,Isolates2,Metagenomics2,cutoff_isolates,cutoff_assigned,cutoff_real)
p1=figure

Isolates=Isolates1;
Metagenomics=Metagenomics1;

AnyIsolatesAtAll=sum(Isolates,1)>0;

Isolates=Isolates(:,AnyIsolatesAtAll);
Metagenomics=Metagenomics(:,AnyIsolatesAtAll);

EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);
%%
Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);
found_either = min([Niso NMetagenomics],[],2)
median_min_found=ceil(median(found_either))
% scatter plot
% Number of lineages in isolates vs metagenomics
subplot(2,1,1)
histogram(found_either,'BinWidth',1,'FaceColor','white')
pbaspect([1 1 1])
xticks([.5:1:10.5])
xticklabels([0:1:10])
xlim([0 10])
ylim([0 5])
hold on
plot([median_min_found+.5 median_min_found+.5],[0 5],'LineWidth',1,'Color','black') % add .5 to put it in the middle of the bin
%%
Isolates=Isolates2;
Metagenomics=Metagenomics2;

AnyIsolatesAtAll=sum(Isolates,1)>0;
MG_excluded_lineages=Metagenomics(:,~AnyIsolatesAtAll);

Isolates=Isolates(:,AnyIsolatesAtAll);
Metagenomics=Metagenomics(:,AnyIsolatesAtAll);

EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);
%%
Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);

% scatter plot
N = numel(Niso);
jitter = rand(N,2)/3;
% Number of lineages in isolates vs metagenomics
Xm=10;
subplot(2,1,2)
% scatter(Niso+jitter(:,1),NMetagenomics+jitter(:,2),'filled','MarkerFaceColor',[0 0 0]);

found_either = min([Niso NMetagenomics],[],2)
median_min_found=median(found_either)
histogram(found_either,'BinWidth',1,'FaceColor','white')
pbaspect([1 1 1])
xticks([.5:1:10.5])
xticklabels([0:1:10])
xlim([0 10])
ylim([0 5])
hold on
plot([median_min_found+.5 median_min_found+.5],[0 5],'LineWidth',1,'Color','black') % add .5 to put it in the middle of the bin
%%



end

