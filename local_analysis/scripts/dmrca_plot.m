%% for plotting scatters
function dmrca_plot(XX,YY,MolecularClockRate,clrsp,str,f)

%% set current figure
set(0,'CurrentFigure',f)

%% set plot limits
nearest_multiple=5;
ylims=[0 nearest_multiple*ceil(max(YY)/nearest_multiple)]; % round up to nearest 5

xlims=[floor(min(XX)) ceil(max(XX))]; 


%% Age vs dMRCA

scatter(XX,YY,'filled',CData=clrsp)
hold on
xlim(xlims)
ylim(ylims)

%% correlation coefficient

[rho,pval]=corr(XX,YY,'type','Pearson');
text(xlims(1),ylims(2)*.95,['R^2=' char(num2str(rho^2))])
text(xlims(1),ylims(2)*.8,['p=' char(num2str(pval))])

%% fit to linear model and plot
p = polyfit(XX, YY, 1);
px = [xlims];
py = polyval(p, px);
plot(px, py, 'LineWidth', 1,'Color','black','LineStyle','--');
xlabel('Age')
ylabel(['dMRCA' make_underscore(str)])

if MolecularClockRate>0
    h=gca;
    yyaxis right
    h.YAxis(2).Limits=h.YAxis(1).Limits/MolecularClockRate;
    h.YAxis(2).TickValues=(0:1:ceil(max(h.YAxis(1).Limits/MolecularClockRate)));
    ylabel(['tMRCA' make_underscore(str)],'Rotation',-90,'VerticalAlignment','middle');

end

end
