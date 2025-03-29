function p=plot_sepi_cacnes_abundance(M,MsubjectTime,clrsgroups)
%% one color/dot for each subject
Cutotypes = M.subject_cutotype;
SIDs = M.SID;
S = numel(SIDs);

clrs = clrsgroups(Cutotypes,:);

%% Get relative species abundance from bracken

SpeciesAbundance=vertcat(MsubjectTime.BrackenAbundanceClustering);
SpeciesNames=MsubjectTime.BrackenNamesClustering(1,:);

CacnesAbundance = SpeciesAbundance(:,SpeciesNames=="Cutibacterium acnes");
SepiAbundance = SpeciesAbundance(:,SpeciesNames=="Staphylococcus epidermidis");

% average across timepoints
CacnesAbundanceSubject = arrayfun(@(x) mean(CacnesAbundance(MsubjectTime.SID==x)), SIDs);
SepiAbundanceSubject = arrayfun(@(x) mean(SepiAbundance(MsubjectTime.SID==x)), SIDs);


%% make plots for both species

p=figure
subplot(2,3,1)
scatter(CacnesAbundanceSubject,M.CacnesOverHuman,'filled','CData',clrs)
ylim([0 1.25])
xlim([0 1])
yticks([0:.25:1.25])
xlabel('C.acnes Rel. Abundance')
ylabel('C.acnes/Human')

subplot(2,3,2)
run_swarmchart(Cutotypes,CacnesAbundanceSubject,clrs,1)
ylabel("C. acnes Rel. abundance")
yticks([.25:.25:1])

subplot(2,3,3)
run_swarmchart(Cutotypes,M.CacnesOverHuman,clrs,1.25)
ylabel("C. acnes/Human")
yticks([0:.25:1.25])

subplot(2,3,4)
scatter(SepiAbundanceSubject,M.SepiOverHuman,'filled','CData',clrs)
xlim([0 .4])
ylim([0 .1])
yticks(.0:.025:.1)
xlabel('S. epi. Rel Abundance')
ylabel('S. epi./Human')
%

subplot(2,3,5)
run_swarmchart(Cutotypes,SepiAbundanceSubject,clrs,.4)
yticks(.0:.1:.4)
ylabel("S. epi. Relative abundance")

subplot(2,3,6)
run_swarmchart(Cutotypes,M.SepiOverHuman,clrs,.1)
yticks([0:.025:.1])
