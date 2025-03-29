function plot_intervals_oh_et_al(AllIntervalData)

% plot all of the types  x of locations seperately, but group p-values all together for FDR
%% Sebaceous sites first
SebaceousLocations=["Al-R" "Ba-R" "Ch-R" "Ea-R" "Mb-R" "Oc-R" "Ra-R"];
[~,SebaceousSitesP_vals_interval, SebaceousSitesP_vals_BrayCurtis]=plot_intervals(AllIntervalData(ismember(AllIntervalData.Location,SebaceousLocations),:));

%% Dry sites
DryLocations=["Hp-R" "Vf-R"];
[~,DrySitesP_vals_interval, DrySitesP_vals_BrayCurtis]=plot_intervals(AllIntervalData(ismember(AllIntervalData.Location,DryLocations),:));
%% Moist locations
MoistLocations=["Ac-R" "Ic-R" "Id-R" "Pc-R"];
[~,MoistSitesP_vals_interval, MoistSitesP_vals_BrayCurtis]=plot_intervals(AllIntervalData(ismember(AllIntervalData.Location,MoistLocations),:));

%% Foot locations
FootLocations=["Ph-R" "Tn-R" "Tw-R"];
[~,FootP_vals_interval, FootP_vals_BrayCurtis]=plot_intervals(AllIntervalData(ismember(AllIntervalData.Location,FootLocations),:));
%%

SitesAll = [SebaceousLocations DryLocations MoistLocations FootLocations];
PVs_BrayCurtis = [SebaceousSitesP_vals_BrayCurtis; DrySitesP_vals_BrayCurtis; MoistSitesP_vals_BrayCurtis; FootP_vals_BrayCurtis];
PVs_Intervals = [SebaceousSitesP_vals_interval; DrySitesP_vals_interval; MoistSitesP_vals_interval; FootP_vals_interval];
SitesWithInsignificantIntervalDifferences = SitesAll(~bh_fdr(PVs_Intervals,.05));
SitesWithInsignificantBrayCurtisDifferences = SitesAll(~bh_fdr(PVs_BrayCurtis,.05));

%%

plot_intervals(AllIntervalData)

end

function [F,PV_int,PV_BC]=plot_intervals(I)

F = figure;
UniqueLocations = unique(I.Location);
L=numel(UniqueLocations);
tiledlayout(L,2)

%%
[PV_int,PV_BC]=deal(zeros(L,1));
subplotidx=1;
for l=1:L
    Il = I(I.Location==UniqueLocations(l),:);

    % first, plot average abundance lost across interval
    nexttile%subplot(L,2,subplotidx);subplotidx=subplotidx+1;
    AbundanceLostCacnes = cellfun(@nanmean,Il.AbundanceLostCacnes);
    AbundanceLostSepi = cellfun(@nanmean,Il.AbundanceLostSepi);
    
    swarmchart(ones(numel(AbundanceLostCacnes),1),AbundanceLostCacnes,'filled','MarkerFaceColor','white','MarkerEdgeColor','black')
    hold on
    swarmchart(2*ones(numel(AbundanceLostSepi),1),AbundanceLostSepi,'filled','MarkerFaceColor','black','MarkerEdgeColor','black')
    xticks([1:2])
    xticklabels(["C. acnes", "S. epidermidis"])
    ylabel('Abundance lost ')
    title(['site:' char(UniqueLocations(l))])
    ylim([0 1])
    yticks([0:.25:1])
    % get p-value
    pv = ranksum(AbundanceLostCacnes,AbundanceLostSepi);
    PV_int(l)=pv;
    text(1.25,.45,['p=' char(string(pv))],'HorizontalAlignment', 'center','VerticalAlignment', 'top');
    xtickangle(0)
    plot([.75 1.25],[median(AbundanceLostCacnes) median(AbundanceLostCacnes)],'LineWidth',1,'Color','black')
    plot([1.75 2.15],[median(AbundanceLostSepi) median(AbundanceLostSepi)],'LineWidth',1,'Color','black')
    % pbaspect([2 1 1])
    h=gca;
    h.LineWidth=.25;

    % then, plot BC changes acoss interval
    nexttile%subplot(L,2,subplotidx);subplotidx=subplotidx+1;
    MaxBrayCurtisCacnes = cellfun(@max,Il.BrayCurtisCacnes);
    maxBrayCurtisSepi = cellfun(@max,Il.BrayCurtisSepi);
    hold on
    for i = 1:size(Il,1)
        plot(Il.CumulativeTime{i},Il.BrayCurtisCacnes{i},'Color',[.85 .85 .85])
        plot(Il.CumulativeTime{i},Il.BrayCurtisSepi{i},'Color','black');
    end
    ylabel('B-C ')
    xlabel('Time')
    title(['site:' char(UniqueLocations(l))])
    ylim([0 1])
    yticks([0:.25:1])
    h=gca;
    h.LineWidth=.25;
    % get p-value
    pv = ranksum(MaxBrayCurtisCacnes,maxBrayCurtisSepi);
    PV_BC(l)=pv;
    if max(horzcat(Il.CumulativeTime{:}))>2
    xlim([0 2.5])
    xticks([0:.5:2.5])
    text(1.25,.45,['p=' char(string(pv))],'HorizontalAlignment', 'center','VerticalAlignment', 'top');
    xtickangle(0)
    else
    xlim([0 .2])
    xticks([0:.05:.2])
    text(.15,.45,['p=' char(string(pv))],'HorizontalAlignment', 'center','VerticalAlignment', 'top');
    xtickangle(0)
    end
    % pbaspect([2 1 1])
end
%%
end