function run_swarmchart(Cutotypes,Y,clrs,Ymax)
%%

%%
swarmchart(Cutotypes,Y,'filled','CData',clrs,'XJitter','density','XJitterWidth',.25)
hold on
arrayfun(@(x) plot([x-.125 x+.125],[median(Y(Cutotypes==x)) median(Y(Cutotypes==x))],'Color','black'), 1:3)
xticklabels(["FC1-children" "FC2-children" "parents"])
%
xlim([.5 3.5])
ylim([0 Ymax])

p12 = char(string(round(ranksum(Y(Cutotypes==1),Y(Cutotypes==2)),4,'significant')));
p23 = char(string(round(ranksum(Y(Cutotypes==2),Y(Cutotypes==3)),4,'significant')));
p13 = char(string(round(ranksum(Y(Cutotypes==1),Y(Cutotypes==3)),4,'significant')));
text(.5,Ymax,['p12=' p12])
text(2.5,Ymax,['p23=' p23])
text(1.5,Ymax*1.1,['p13=' p13])
%