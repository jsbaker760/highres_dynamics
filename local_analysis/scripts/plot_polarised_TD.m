function LD=plot_polarised_TD(LD,MinNiso)

LD.dmrca_vs_pairwise_d = arrayfun(@(x) PolarizedTajimaD_lineage(LD.has_mutation_ingroup_goodpos_no_reco{x},MinNiso),1:size(LD,1))';

pv=ranksum(LD.dmrca_vs_pairwise_d(LD.SpeciesName=="cacnes"),LD.dmrca_vs_pairwise_d(LD.SpeciesName~="cacnes"));

f=figure;
subplot(1,2,1);
histogram(LD.dmrca_vs_pairwise_d(LD.SpeciesName=="cacnes"),'BinWidth',.1,'Normalization','probability','FaceColor','white');
title("C. acnes")
xlabel('µ pairwise distance ÷ 2 µdMRCA')
ylabel('Proportion lineages')
pbaspect([1 1 1]);xticks([0 .5 1])
ylim([0 .3]);yticks([0 .1 .2 .3])
subplot(1,2,2);
histogram(LD.dmrca_vs_pairwise_d(LD.SpeciesName=="sepi"),'BinWidth',.1,'Normalization','probability','FaceColor','black');
title("S. epidermidis")
ylim([0 .3]);yticks([0 .1 .2 .3])
pbaspect([1 1 1]);xticks([0 .5 1])
xlabel('µ pairwise distance ÷ 2 µdMRCA')
ylabel('Proportion lineages')