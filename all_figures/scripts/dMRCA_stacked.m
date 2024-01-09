function p= dMRCA_stacked(X,Y,Cutotypes,clrs,MolecularClockRate,ttl)

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
[N1, edges1, index1]=histcounts(Y(B1),0:step:MaxBin);
[N2, edges2, index2]=histcounts(Y(B2),0:step:MaxBin);
[N3, edges3, index3]=histcounts(Y(B3),0:step:MaxBin);

Nsteps=numel(N1);
Nall = [N1;N2;N3]
Cdata = [.337 .705 .913;0 .4471 .698;.8253 .3686 .0039];

b=bar(Nall','stacked','FaceColor','flat');
for i = 1:3
    b(i).CData=Cdata(i,:);
end
% xlim([0 Nsteps+1]);
xticks([1:Nsteps])
xticklabels([edges1(2:end)])
xlabel('dMRCA')
ylabel('Number of phylogenies')
if MolecularClockRate>0
    hAx(1)=gca;
    hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','yAxisLocation','right', 'color','none');
    hold(hAx(2),'on')
    b=bar(Nall','stacked','FaceColor','flat');
    for i = 1:3
        b(i).CData=Cdata(i,:);
    end
    foo=1;
    % xlim([.5 Nsteps+.5]);
    xticks([1:Nsteps])
    xticklabels(round([edges1(2:end)./MolecularClockRate],1,'decimals'))
    hold on
    xlabel('tMRCA')
    ylabel('Number of phylogenies')
end
%%
% if MolecularClockRate>0
% 
%     % set(0,'CurrentFigure')
%     Y=Y./MolecularClockRate;
%     %%
%     subplot(2,1,2)
%     step=1;
%     MaxBin=max(Y);
%     MaxBin = MaxBin-mod(MaxBin,step)+step;
% 
%     B1=Cutotypes==1;
%     B2=Cutotypes==2;
%     B3=Cutotypes==3;
%     [N1, edges1, index1]=histcounts(Y(B1),0:step:MaxBin);
%     [N2, edges2, index2]=histcounts(Y(B2),0:step:MaxBin);
%     [N3, edges3, index3]=histcounts(Y(B3),0:step:MaxBin);
% 
%     Nsteps=numel(N1);
%     Nall = [N1;N2;N3]
%     Cdata = [.337 .705 .913;0 .4471 .698;.8253 .3686 .0039];
% 
%     b=bar(Nall','stacked','FaceColor','flat');
%     for i = 1:3
%         b(i).CData=Cdata(i,:);
%     end
%     xlim([0 Nsteps+1]);
%     xticks([1:Nsteps])
%     xticklabels([edges1(2:end)])
%     hold on
%     xlabel('tMRCA')
%     ylabel('Number of phylogenies')
% end