function p=plot_sepi_cacnes_abundance(M,clrs_groups)

%%

Cutotypes = M.subject_cutotype;

%% one color/dot for each subject

clrs = clrs_groups(Cutotypes,:);

%% Get relative speciesabundance from bracken

SpeciesAbundance=vertcat(M.BrackenSpecies{:});
CacnesAbundance = SpeciesAbundance(:,M.NamesSpecies(1,:)=="Cutibacterium acnes");
SepiAbundance = SpeciesAbundance(:,M.NamesSpecies(1,:)=="Staphylococcus epidermidis");

%% make plots for both species

p=figure
subplot(2,3,1)
scatter(CacnesAbundance,M.CacnesOverHuman,'filled','CData',clrs)
ylim([0 1.25])
xlim([0 1])
yticks([0:.25:1.25])
xlabel('C.acnes Rel. Abundance')
ylabel('C.acnes/Human')

subplot(2,3,2)
run_swarmchart(Cutotypes,CacnesAbundance,clrs,1)
ylabel("C. acnes Rel. abundance")
yticks([.25:.25:1])

subplot(2,3,3)
run_swarmchart(Cutotypes,M.CacnesOverHuman,clrs,1.25)
ylabel("C. acnes/Human")
yticks([0:.25:1.25])

subplot(2,3,4)
scatter(SepiAbundance,M.SepiOverHuman,'filled','CData',clrs)
xlim([0 .4])
ylim([0 .1])
yticks(.0:.025:.1)
xlabel('S. epi. Rel Abundance')
ylabel('S. epi./Human')
%

subplot(2,3,5)
run_swarmchart(Cutotypes,SepiAbundance,clrs,.4)
yticks(.0:.1:.4)
ylabel("S. epi. Relative abundance")

subplot(2,3,6)
run_swarmchart(Cutotypes,M.SepiOverHuman,clrs,.1)
yticks([0:.025:.1])

%%