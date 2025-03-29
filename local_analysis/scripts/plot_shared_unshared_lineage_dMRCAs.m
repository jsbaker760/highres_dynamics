function plot_shared_unshared_lineage_dMRCAs(LineageDMRCAs,ManuscriptColors)

%% split first for C. acnes
figure
subplot(1,2,1)
L = LineageDMRCAs(LineageDMRCAs.Nisolates>=5&LineageDMRCAs.Species=="cacnes",:);
average_dMRCA = cellfun(@mean ,L.dMRCA_subject_isolates);
unshared = L.Unshared;LUnshared=L(unshared,:);LShared = L(~unshared,:);
% plot difference between shared and unshared
swarmchart(unshared+1, average_dMRCA,'filled','CData',ManuscriptColors.CutotypeColors(L.Cutotype,:),'XJitterWidth',.5)
hold on
plot([.75 1.25],[median(average_dMRCA(unshared==0)) median(average_dMRCA(unshared==0))],'LineWidth',2,'Color','black')
plot([1.75 2.25],[median(average_dMRCA(unshared==1)) median(average_dMRCA(unshared==1))],'LineWidth',2,'Color','black')
xticks([1:2])
xticklabels(["Shared" "Unshared"])
ylabel('average dMRCA')
title('C. acnes')
pbaspect([1 1 1])
p=ranksum(average_dMRCA(unshared),average_dMRCA(~unshared));
text(1,30,['p=' char(string(p))])
%
L = LineageDMRCAs(LineageDMRCAs.Nisolates>=5&LineageDMRCAs.Species=="sepi",:);
subplot(1,2,2)
clrsgroups=ManuscriptColors.CutotypeColors;
average_dMRCA = cellfun(@mean ,L.dMRCA_subject_isolates);
unshared = L.Unshared;
LUnshared=L(unshared,:);
LShared = L(~unshared,:);
swarmchart(unshared+1, average_dMRCA,'filled','CData',ManuscriptColors.CutotypeColors(L.Cutotype,:),'XJitterWidth',.5)
hold on
plot([.75 1.25],[median(average_dMRCA(unshared==0)) median(average_dMRCA(unshared==0))],'LineWidth',2,'Color','black')
plot([1.75 2.25],[median(average_dMRCA(unshared==1)) median(average_dMRCA(unshared==1))],'LineWidth',2,'Color','black')

xticks([1:2])
xticklabels(["Shared" "Unshared"])
ylabel('average dMRCA')
p=ranksum(average_dMRCA(unshared),average_dMRCA(~unshared));
text(1,30,['p=' char(string(p))])
pbaspect([1 1 1])
title('S. epidermidis')