function p= dMRCA_histograms(X,Y,Cutotypes)

p = figure();
% if MolecularClockRate>0
%     subplot(2,1,1)
% end
% set(0,'CurrentFigure')

step=5;
MaxBin=max(Y);
MaxBin = MaxBin-mod(MaxBin,step)+step;

B1=Cutotypes==1;
B2=Cutotypes==2;
B3=Cutotypes==3;
[N1, ~, ~]=histcounts(Y(B1),0:step:MaxBin);
[N2, ~, ~]=histcounts(Y(B2),0:step:MaxBin);
[N3, ~, ~]=histcounts(Y(B3),0:step:MaxBin);

Nsteps=numel(N1);
Nall = [N1;N2;N3];
Cdata = [.337 .705 .913;0 .4471 .698;.8253 .3686 .0039];

for i = 1:3
    subplot(3,1,i)
    D = Y(Cutotypes==i);
    histogram(D,'FaceColor',Cdata(i,:),'Normalization','probability','BinWidth',step);
    hold on
    xlim([0 MaxBin])
    xticklabels([0:5:MaxBin])
    xlabel('dMRCA')
    ylabel('Proportion of phylogenies')
    ylim([0 .8])

    yticks([0:.2:.8])
    plot([median(D) median(D)],[0 1],'LineWidth',1,'Color','black')
    pbaspect([1 1 1])
    hAx(1)=gca;
end

%%
p_all = ranksum(Y(Cutotypes==3),Y(Cutotypes==1|Cutotypes==2));

['P=' char(string(p_all)) '(parents vs children)']

