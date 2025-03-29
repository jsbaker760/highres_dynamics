%% To plot everything
function p= dMRCA_regression(X,Yi,Cutotypes,clrs,MolecularClockRate,str)

p = figure;
% plot for everybody
hold on
dmrca_plot(X,Yi,MolecularClockRate,clrs,str,p)
pbaspect([2 1 1])
set(0,'CurrentFigure')
% plot for children
subplot(2,3,1)
hold on
XX=X(Cutotypes<3);
YY=Yi(Cutotypes<3);
clrsp=clrs(Cutotypes<3,:);
dmrca_plot(XX,YY,MolecularClockRate,clrsp,str,p)
pbaspect([1 1 1])
% plot for parents
subplot(2,3,2)
hold on
XX=X(Cutotypes==3);
YY=Yi(Cutotypes==3);
clrsp=clrs(Cutotypes==3,:);
dmrca_plot(XX,YY,MolecularClockRate,clrsp,str,p)
pbaspect([1 1 1])
% plot for everybody
subplot(2,3,4:5)
hold on
dmrca_plot(X,Yi,MolecularClockRate,clrs,str,p)
pbaspect([2 1 1])

% swamrmcharts
subplot(2,3,3)
X = Cutotypes;
lbls=["FC1-Children" "FC2-Children","Parents"];
swarms(X,Yi,clrs,lbls,MolecularClockRate,str,p)
subplot(2,3,6)
X = Cutotypes;
X(X==2)=1;X(X==3)=2;
lbls=["All Children" "Parents"];
swarms(X,Yi,clrs,lbls,MolecularClockRate,str,p)

end