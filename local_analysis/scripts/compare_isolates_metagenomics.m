function [CacnesLineagesEitherDataSource,SepiLineagesEitherDataSource]=compare_isolates_metagenomics(IT,M,cutoff_isolates,cutoff_assigned,cutoff_real)
%% Get abundances to use in comparison for both species
plothandles = struct;
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

%% plot coverage difference vs recall
coverage=M.SepiMGCoverage;
[plothandles.sepi_coverage_dependence]= plot_coverage_difference(SepiIsolateNumber,MetagenomicsSepiAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.25 .25 .25],coverage);
ax=gca;
ax.Title.String='S. epidermidis';
coverage=M.CacnesMGCoverage;
[plothandles.cacnes_coverage_dependence]= plot_coverage_difference(CacnesIsolateNumber,MetagenomicsCacnesAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.25 .25 .25],coverage);
ax=gca;
ax.Title.String='C. acnes';
%% Corrospondence between MG and isolates
[plothandles.CacnesComparisons]= plot_comparisons(CacnesIsolateNumber,MetagenomicsCacnesAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.85 .85 .85]);
ax=gca;
ax.Parent.Children(3).Title.String='C. acnes';
[plothandles.SepiComparisons]= plot_comparisons(SepiIsolateNumber,MetagenomicsSepiAbundance,cutoff_isolates,cutoff_assigned,cutoff_real,[.25 .25 .25]);
ax=gca;
ax.Parent.Children(3).Title.String='S. epidermidis';

%% Find isolates found in EITHER data source (union)

CacnesLineagesEitherDataSource = CacnesIsolateNumber>0|MetagenomicsCacnesAbundance>cutoff_real;
CacnesLineagesEitherDataSource(sum(CacnesIsolateNumber,2)<cutoff_isolates|sum(MetagenomicsCacnesAbundance,2)<cutoff_assigned,:)=0;
SepiLineagesEitherDataSource = SepiIsolateNumber>0|MetagenomicsSepiAbundance>cutoff_real;
SepiLineagesEitherDataSource(sum(SepiIsolateNumber,2)<cutoff_isolates|sum(MetagenomicsSepiAbundance,2)<cutoff_assigned,:)=0;

f = figure;
subplot(2,1,1)
NCacnesEither = sum(CacnesLineagesEitherDataSource,2);
NCacnesEither(NCacnesEither==0)=[];
histogram(NCacnesEither,'BinWidth',.99,'FaceColor','white')
hold on
plot([median(NCacnesEither)+.5 median(NCacnesEither)+.5],[0 10],'LineWidth',2,'Color','black')
xticks(.5:2:12.5)
xticklabels(0:2:12)
xlim([0.5 12.5])
ylim([0 6])
xlabel('Number of coexisting lineages')
ylabel('Number of samples')
title('C. acnes')
subplot(2,1,2)
NSepiEither = sum(SepiLineagesEitherDataSource,2);
NSepiEither(NSepiEither==0)=[];
histogram(NSepiEither,'BinWidth',.99,'FaceColor','black')
hold on
plot([median(NSepiEither)+.5 median(NSepiEither)+.5],[0 10],'LineWidth',2,'Color','black')
xticks(.5:2:12.5)
xticklabels(0:2:12)
xlim([0.5 12.5])
ylim([0 6])
xlabel('Number of coexisting lineages')
ylabel('Number of samples')
title('S. epidermidis')
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


function [p]= plot_comparisons(Isolates,Metagenomics,cutoff_isolates,cutoff_assigned,cutoff_real,clr)
p=figure;

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



%% make cutoffs and jitter

Metagenomics(Metagenomics<cutoff_real)=0;

% Number of clades found in both
Niso=sum(Isolates>0,2);
NMetagenomics=sum(Metagenomics>0,2);

% scatter plot
N = numel(Niso);
jitter = rand(N,2)/3;

%% Number of lineages in isolates vs metagenomics
Xm=10;
subplot(3,1,1)
scatter(Niso+jitter(:,1),NMetagenomics+jitter(:,2),'filled','MarkerFaceColor',clr);
hold on
plot([0 Xm],[0 Xm],'Color','black')
xlim([0 Xm])
ylim([0 Xm])
pbaspect([1 1 1])
[rho,p]=corr(Niso,NMetagenomics,'type','Pearson');
text(0,Xm,['p=' char(string(p))])
text(0,Xm-2,['R2=' char(string(rho^2))])

xlabel('Number of lineages (Isolates)')
ylabel('Number of lineages (Metagenomics)')

%% Abundance of lineages found at both tps

FoundBoth=(Isolates>0&Metagenomics>0);
MetagenomicsBoth=Metagenomics(FoundBoth);
IsolatesBoth=Isolates(FoundBoth);
subplot(3,1,2)

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

subplot(3,1,3)

MetagenomicsMissingInIso=Metagenomics(Metagenomics>0&Isolates==0);
histogram(MetagenomicsMissingInIso,'BinWidth',.01)
xlim([0 1])
% ylim([0 1])
pbaspect([1 1 1])
xlabel('Relative abundance (Metagenomics) of lineages not found in isolates and found in metagenomics')
ylabel('number of instances')
xticks([0:.2:1])

%%




end


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
found_either = min([Niso NMetagenomics],[],2);
median_min_found=ceil(median(found_either));
% scatter plot
% Number of lineages in isolates vs metagenomics
%%
found_either_cacnes = found_either;
%%
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
found_either_sepi = found_either;
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



%%

function [p]= plot_coverage_difference(Isolates,Metagenomics,cutoff_isolates,cutoff_assigned,cutoff_real,clr,coverage)

%% remove rows with insufficient data


EnoughIsolates = sum(Isolates,2)>cutoff_isolates;
EnoughMG = sum(Metagenomics,2)>cutoff_assigned;

EnoughBoth=EnoughIsolates&EnoughMG;

Isolates=Isolates(EnoughBoth,:);
Metagenomics=Metagenomics(EnoughBoth,:);
coverage=coverage(EnoughBoth,:);

% normalize isolates
Isolates=Isolates./sum(Isolates,2);

%%

Metagenomics(Metagenomics<cutoff_real)=0;

IsoMissingInMetagenomics=Isolates>0&Metagenomics==0;
% foo=1;

FoundEither=Isolates>0|Metagenomics>0;

% Y = sum(IsoMissingInMetagenomics,2)./sum(FoundEither>0,2);

IsoMissingInMetagenomics = sum(IsoMissingInMetagenomics,2);
Isolates=sum(Isolates>0,2);

Y = IsoMissingInMetagenomics./Isolates;
p=figure;

%%

[f,gof] = fit(coverage,Y,'exp1');
plot(f)
text(20,.8,['adjR2=' string(gof.adjrsquare)])
[~,pv]=corr(coverage,Y,'Type','Spearman');
text(20, .9,string(pv));
hold on
scatter(coverage,Y,'filled')
ylim([0 1])
xlim([0 100])
yticks([0:.2:1])
xlabel('proportion of isolate lineages not found in metagenomes')
ylabel('Reference coverage')
end