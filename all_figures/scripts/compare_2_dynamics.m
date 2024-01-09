function p=compare_2_dynamics(D1,D2,type)

%%
load data/ColorBlind.mat
%%
[Nlost1,Nextant1,Alost1]=get_lost(D1)
[Nlost2,Nextant2,Alost2]=get_lost(D2)

%%

Plost1 = arrayfun(@(x) mean(Nlost1{x}'./Nextant1{x}), 1:numel(Nlost1));
Plost2 = arrayfun(@(x) mean(Nlost2{x}'./Nextant2{x}), 1:numel(Nlost2));

if size(Plost2,1)>1
    Plost2=Plost2'
end
if size(Plost1,1)>1
    Plost1=Plost1'
end
%%
clrs3 = [ColorBlind.LightBlue ; ColorBlind.Blue ; ColorBlind.Vermillion];

clrs = [cellfun(@unique,D1.CutoTPs) ; cellfun(@unique,D2.CutoTPs)];
clrs = clrs3(clrs,:)
%%
p=figure
subplot(1,2,1)
swarmchart([ones(numel(Plost1),1); 2*ones(numel(Plost2),1)],[Plost1 Plost2],'filled','XJitter','density','XJitterWidth',.5,'CData',clrs)
xticks(1:2)
xticklabels(["C. acnes" "S. epi"])
ylabel(strjoin(["Proportion" type "lost (average)"]))
ylim([0 1])
text(1.5,1,string(round(ranksum(Plost1,Plost2),2,'significant')))
pbaspect([1 1 1])
hold on
plot([.75 1.25], [median(Plost1) median(Plost1)],'Color','black')
plot([1.75 2.25], [median(Plost2) median(Plost2)],'Color','black')
subplot(1,2,2)
swarmchart([ones(numel(Plost1),1); 2*ones(numel(Plost2),1)],[cellfun(@mean,Alost1); cellfun(@mean,Alost2)],'filled','XJitter','density','XJitterWidth',.5,'CData',clrs)
xticks([1:3])
xticklabels(["C. acnes" "S. epi"])
ylim([0 1])
ylabel(strjoin(["Abundance" type "lost (average)"]))
text(1.5,1,string(round(ranksum(cellfun(@mean,Alost1),cellfun(@mean,Alost2)),2,'significant')))
pbaspect([1 1 1])
hold on
plot([.75 1.25], [median(cellfun(@mean,Alost1)) median(cellfun(@mean,Alost1))],'Color','black')
plot([1.75 2.25], [median(cellfun(@mean,Alost2)) median(cellfun(@mean,Alost2))],'Color','black')

%%
function [Nlost,Nextant,Alost]=get_lost(D)

R = size(D,1);
[Nlost,Nextant,Alost]=deal(cell(R,1));
for i = 1:size(D,1)
    A = D.Abundances{i};
    N=sum(A>0,2);
    Nextant{i}=N(1:end-1);
    L = arrayfun(@(x) sum(A(x,:)>0&(A(x+1,:)==0)), 1:(size(A,1)-1));
    Nlost{i}=L;
    Al = arrayfun(@(x) sum(A(x,A(x,:)>0&(A(x+1,:)==0))), 1:(size(A,1)-1));
    Alost{i}=mean(Al);
end