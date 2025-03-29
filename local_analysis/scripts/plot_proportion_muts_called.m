function plot_proportion_muts_called(LineageData)
%% Figure S24
LineageData_plot=LineageData;

% remove lineages without de novo SNPs 
em=cellfun(@isempty,LineageData.goodpos);
LineageData_plot(em,:)=[];
goodpos_ingroup_calls = arrayfun(@(x) {LineageData_plot.calls_analysis{x}(:,~LineageData_plot.outgroup{x})>0} , 1:size(LineageData_plot,1));
proportion_calls = cellfun(@(x) {mean(x,1)}, goodpos_ingroup_calls);

% plot
cacnes_calls = horzcat(proportion_calls{LineageData_plot.SpeciesName=="cacnes"});
sepi_calls = horzcat(proportion_calls{LineageData_plot.SpeciesName=="sepi"});
f=figure;
subplot(2,1,1)
histogram(cacnes_calls,'BinWidth',.01,'FaceColor','white')
xlabel("Proportion within-lineage mutation positions called")
ylabel("number of isolates")
xlim([.5 1])
ylim([0 2000])
subplot(2,1,2)
histogram(sepi_calls,'BinWidth',.01,'FaceColor','black')
xlabel("Proportion within-lineage mutation positions called")
ylabel("number of isolates")
xlim([.5 1])
ylim([0 2000])
%