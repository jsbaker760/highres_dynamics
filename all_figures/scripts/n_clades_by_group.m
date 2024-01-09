function [p,N,Tbl] = n_clades_by_group(Tbl,A,min_assigned,cutoff_real,ttl)

%% Get the amount assigned

Assigned = sum(A,2);
enough=Assigned>min_assigned;

%% pare down data to rows with enough assigned

Tbl=Tbl(enough,:);
A=A(enough,:);
A(A<cutoff_real)=0;

%% Number of clades

N = sum(A>0,2);

%% get average N per subject, as well as subject's cutotype

[uSID,ia,~]=unique(Tbl.SID);
Cutotypes = Tbl.subject_cutotype(ia);
Cutotypes(contains(uSID,"P"))=3;
Nmean = arrayfun(@(x) mean(N(Tbl.SID==x)), uSID);

%% get colors for plot

load data/ColorBlind.mat
clrs=[ColorBlind.LightBlue; ColorBlind.Blue; ColorBlind.Vermillion];
clrs=clrs(Cutotypes,:);

%% group daya

DsGrouped = arrayfun(@(x) {Nmean(Cutotypes==x)}, 1:3);
clrsGroups= arrayfun(@(x) {clrs(Cutotypes==x,:)}, 1:3);

AllChildren = find(Cutotypes<3);
DsAllChildren = Nmean(AllChildren,:);
clrsAllChildren = clrs(AllChildren,:);

Ds=[DsGrouped(1) DsGrouped(2) DsAllChildren DsGrouped(3)];
clrs = [clrsGroups(1) clrsGroups(2) clrsAllChildren clrsGroups(3)];

Xs = arrayfun(@(x) {repmat(x,numel(Ds{x}),1)}, 1:4);
Xs=vertcat(Xs{:});
Ys=vertcat(Ds{:});
clrs=vertcat(clrs{:});
swarmchart(Xs,Ys,'filled','CData',clrs,'XJitter','rand','XJitterWidth',.25);
hold on
arrayfun(@(x) plot([x-.2 x+.2],[median(Ds{x}) median(Ds{x})],'LineWidth',1,'Color','black'), [1:4])

%% format

Ymax = 12%ceil(max(Nmean));
ylim([0 Ymax]);
xlim([.5 4.5]);
xticks(1:4)
xticklabels(["FC1-Children" "FC2-Children" "AllChildren","Parents"])
title(ttl)
%% calculate p-values

p12=char(string(round(ranksum(Ds{1},Ds{2}),2,'significant')));
p24=char(string(round(ranksum(Ds{2},Ds{4}),2,'significant')));
p14=char(string(round(ranksum(Ds{1},Ds{4}),2,'significant')));
pPvsC=char(string(round(ranksum(Ds{3},Ds{4}),2,'significant')));

% [~,p12,~]=kstest2(Ds{1},Ds{2})
% [~,p24,~]=kstest2(Ds{4},Ds{2})
% [~,p14,~]=kstest2(Ds{1},Ds{4})
% [~,pPvsC,~]=kstest2(Ds{3},Ds{4})
% 
% p12=char(string(round(p12,2,'significant')));
% p24=char(string(round(p24,2,'significant')));
% p14=char(string(round(p14,2,'significant')));
% pPvsC=char(string(round(pPvsC,2,'significant')));

%% add p-values
text(.5,Ymax,['p12=' p12])
text(.5,Ymax-1,['p23=' p24])
text(.5,Ymax-2,['p13=' p14])
text(.5,Ymax-3,['pParentsChildren=' pPvsC])
%%

p=gca;

%%

scatter(Tbl.Age(ia),Nmean)
