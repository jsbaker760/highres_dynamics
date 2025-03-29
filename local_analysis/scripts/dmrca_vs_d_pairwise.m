function [LineageDMRCAs] = dmrca_vs_d_pairwise(LineageDMRCAs,ManuscriptColors)

% split into two structures, one for C. acnes and one for S. epidermidis to
% run function seperately on each
LCacnes = LineageDMRCAs(LineageDMRCAs.Species=="cacnes",:);
LSepi = LineageDMRCAs(LineageDMRCAs.Species=="sepi",:);

clrs3=ManuscriptColors.CutotypeColors;

f=figure
subplot(2,2,1)
[CacnesTajimaD,CacnesN]=PolarizedTajimaD_subjects(LCacnes.SubjectTreeSNPs_isolates,10);
pbaspect([1 1 1])
title('C. acnes')

subplot(2,2,3)
[SepiTajimaD,SepiN]=PolarizedTajimaD_subjects(LSepi.SubjectTreeSNPs_isolates,10);
pbaspect([1 1 1])
title('S. epidermidis')
subplot(2,2,2)
pbaspect([1 1 1])
scatter(CacnesN,CacnesTajimaD,'filled','CData',clrs3(LCacnes.Cutotype,:))
xlabel('N isolates')
ylabel('Polarized Tajima D')
xlim([0 200]);
title('C. acnes')

subplot(2,2,4)
scatter(SepiN,SepiTajimaD,'filled','CData',clrs3(LSepi.Cutotype,:))
xlabel('N isolates')
ylabel('Polarized Tajima D')
xlim([0 200]);
pbaspect([1 1 1])
title('S. epidermidis')

%%

f2=figure

subplot(2,2,1)
X=LSepi.Nisolates;
Y=cellfun(@mean,LSepi.dMRCA_subject_isolates);
scatter(X,Y,'filled','CData',clrs3(LSepi.Cutotype,:))
hold on
xlim([0 200])
ylim([0 60])
ft1 = fit(X,Y,'poly1')
[R,pv]=corr(X,Y)
eq = ['y=' char(string(round(ft1.p1,2,'significant'))) 'x+' char(string(round(ft1.p2,2,'significant')))]
plot(ft1)
l=legend
l.String={eq}
text(20,58,['R^2=' char(string(round(R^2,2,'significant')))])
text(20,50,['P=' char(string(round(pv,2,'significant')))])
title('S. epidermidis')
xlabel('Number of isolates (subject)')
ylabel('dMRCA (subject)')
%
subplot(2,2,3)
X=LCacnes.Nisolates;
Y=cellfun(@mean,LCacnes.dMRCA_subject_isolates);
scatter(X,Y,'filled','CData',clrs3(LCacnes.Cutotype,:))
hold on
xlim([0 200])
ylim([0 60])
ft1 = fit(X,Y,'poly1')
[R,pv]=corr(X,Y)
eq = ['y=' char(string(round(ft1.p1,2,'significant'))) 'x+' char(string(round(ft1.p2,2,'significant')))]
plot(ft1)
l=legend
l.String={eq}
text(20,58,['R^2=' char(string(round(R^2,2,'significant')))])
text(20,50,['P=' char(string(round(pv,2,'significant')))])
xlabel('Number of isolates (subject)')
ylabel('dMRCA (subject)')
title('C. acnes')

subplot(2,2,2)
X=LSepi.Nisolates;
Y=cellfun(@mean,LSepi.dMRCA_subject_isolates);
scatter(X,Y,'filled','CData',clrs3(LSepi.Cutotype,:))
hold on
xlim([0 40])
ylim([0 30])
ft1 = fit(X,Y,'poly1')
[R,pv]=corr(X,Y)
eq = ['y=' char(string(round(ft1.p1,2,'significant'))) 'x+' char(string(round(ft1.p2,2,'significant')))]
plot(ft1)
l=legend
l.String={eq}
text(20,58,['R^2=' char(string(round(R^2,2,'significant')))])
text(20,50,['P=' char(string(round(pv,2,'significant')))])
title('S. epidermidis')
xlabel('Number of isolates (subject)')
ylabel('dMRCA (subject)')
%
subplot(2,2,4)
X=LCacnes.Nisolates;
Y=cellfun(@mean,LCacnes.dMRCA_subject_isolates);
scatter(X,Y,'filled','CData',clrs3(LCacnes.Cutotype,:))
hold on
xlim([0 40])
ylim([0 30])
ft1 = fit(X,Y,'poly1')
[R,pv]=corr(X,Y)
eq = ['y=' char(string(round(ft1.p1,2,'significant'))) 'x+' char(string(round(ft1.p2,2,'significant')))]
plot(ft1)
l=legend
l.String={eq}
text(20,58,['R^2=' char(string(round(R^2,2,'significant')))])
text(20,50,['P=' char(string(round(pv,2,'significant')))])
xlabel('Number of isolates (subject)')
ylabel('dMRCA (subject)')
title('C. acnes')
