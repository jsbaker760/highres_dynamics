%% for plotting swarms
function p=swarms(X,Y,clrs,lbls,MolecularClockRate,str,f)

set(0,'CurrentFigure',f)
hold on
nearest_multiple=5;
xlims=[0 ceil(max(X))+.5]; 
ylims=[0 nearest_multiple*ceil(max(Y)/nearest_multiple)]; % round up to nearest 5

%% parse data into cells for p values later

Ds = arrayfun(@(x) {Y(X==x)}, 1:max(X));

%% make swamrchart

swarmchart(X,Y,'filled','XJitter','density','CData',clrs,'XJitterWidth',.2)
hold on
% add median lines
arrayfun(@(x) plot([x-.2 x+.2], [median(Y(X==x)) median(Y(X==x))],'LineStyle','-','Color','black','LineWidth',1), double(1:max(X)))

pbaspect([1 1 1])
% add labels
xlim(xlims)
xticks(1:max(X))
xticklabels(lbls);
% y label
ylim(ylims)
ylabel(['dMRCA' make_underscore(str)])


% add p values
if numel(Ds)==2
    p=char(num2str(round(ranksum(Ds{1},Ds{2}),2,"significant")));
    text(.5,ylims(2),['p=' p])
elseif numel(Ds)==3
        p12 = char(num2str(round(ranksum(Ds{1},Ds{2}),2,"significant")));
        text(.5,ylims(2),['p12=' p12])
        p13 = char(num2str(round(ranksum(Ds{1},Ds{3}),2,"significant")));
        text(1,ylims(2)*1.1,['p13=' p13])
        p23 = char(num2str(round(ranksum(Ds{3},Ds{2}),2,"significant")));
        text(2.5,ylims(2),['p23=' p23])
        p4 = char(num2str(round(ranksum(Ds{3},[Ds{2} ; Ds{1}]),2,"significant")));
        text(3,ylims(2)*1.1,['p4=' p4])        
else
    error('too many input Cutotypes')
end

% add tMRCA axis

if MolecularClockRate>0
    h=gca;
    yyaxis right
    h.YAxis(2).Limits=h.YAxis(1).Limits/MolecularClockRate;
    h.YAxis(2).TickValues=(0:1:ceil(max(h.YAxis(1).Limits/MolecularClockRate)));
    ylabel(['tMRCA' make_underscore(str)],'Rotation',-90,'VerticalAlignment','middle');
end

end