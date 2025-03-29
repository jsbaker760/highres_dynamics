function f=alpha_diversity_subject(M,clrs)
%% get abundances

[SID,idx] = unique(M.SID);
A = arrayfun(@(x) {mean(M.BrackenAbundanceClustering(M.SID==x,:),1)}, SID);
A = vertcat(A{:});
A = A./sum(A,2);

%% get colors
Cutotypes = M.subject_cutotype(idx);
clrs = clrs(Cutotypes,:);

%% calculate simpson
AlphaDiv= arrayfun(@(x) calc_alpha_div(A(x,:),'simpson'), 1:size(A,1));

% make plot with swarmchat
jw=.25;

f=figure;
hold on
swarmchart(Cutotypes,AlphaDiv,'filled','CData',clrs,'XJitter','density','XJitterWidth',jw)
arrayfun(@(x) plot([x-jw x+jw],[median(AlphaDiv(Cutotypes==x)) median(AlphaDiv(Cutotypes==x))],'Color','black','LineWidth',1), 1:3)

xlim([.5 3.5])
ylim([0 1])
yticks([0:.25:1])
xticks(1:3)
xticklabels(["FC1-Children","FC2-Children","Parents"])

% Get p-values
p12 = ranksum(AlphaDiv(Cutotypes==1),AlphaDiv(Cutotypes==2));
p23 = ranksum(AlphaDiv(Cutotypes==3),AlphaDiv(Cutotypes==2));
p13 = ranksum(AlphaDiv(Cutotypes==1),AlphaDiv(Cutotypes==3));
%
text(.5,1,['p12=' char(string(round(p12,2,'significant')))])
text(1.5,1,['p23=' char(string(round(p23,2,'significant')))])
text(2.5,1,['p13=' char(string(round(p13,2,'significant')))])

