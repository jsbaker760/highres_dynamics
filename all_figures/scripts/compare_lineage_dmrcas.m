function [fMin,fMean,fMax,fAll,fStacked]=compare_lineage_dmrcas(L,MolecularClockRate,metricname,clrsgroups,str)
% metric is the data to measure, can be by isolate or genotype
% with subtly different but similar results

metric = L.(metricname);

%% turn into array if a cell, and average along the first dimention

if iscell(metric)
    metric = cellfun(@mean,metric);
else
    metric = arrayfun(@mean,metric);
end

%% get unique subjects and their cutotypes

[uSID,uidx] = unique(L.SID);
Cutotypes = L.Cutotype(uidx);
% number of subjects

%% colors for each dot

clrs_subjects=clrsgroups(Cutotypes,:);


%% get the age of subject (X) and the mean, min, max Y values per subject 

X = arrayfun(@(x) mean(L.Age(x==L.SID)), uSID);
Ymin = arrayfun(@(x) min(metric(x==L.SID)), uSID);
Ymean = arrayfun(@(x) mean(metric(x==L.SID)), uSID);
Ymax = arrayfun(@(x) max(metric(x==L.SID)), uSID);

%% plot min, mean, max dMRCA/tMRCA, with one dot per subject

fMin=dMRCA_regression(X,Ymin,Cutotypes,clrs_subjects,MolecularClockRate,['min']);
fMean=dMRCA_regression(X,Ymean,Cutotypes,clrs_subjects,MolecularClockRate,['mean']);
fMax=dMRCA_regression(X,Ymax,Cutotypes,clrs_subjects,MolecularClockRate,['max']);

%%
if MolecularClockRate>0
    TmrcaAll=metric./MolecularClockRate
    max_molclock_rate = max(TmrcaAll)
    median_molclock_rate=median(TmrcaAll)
    tMRCA_by_group_all=arrayfun(@(x) {TmrcaAll(L.Cutotype==x)}, 1:3)
    median_by_group=cellfun(@median,tMRCA_by_group_all)
    tMRCAmin=Ymin./MolecularClockRate;
    tMRCAmin_by_group_all=arrayfun(@(x) {tMRCAmin(Cutotypes==x)}, 1:3)
    median_by_group=cellfun(@median,tMRCAmin_by_group_all)
    mean_by_group=cellfun(@mean,tMRCAmin_by_group_all)

end
%% colors for all dots

clrs_all=clrsgroups(L.Cutotype,:);

%% as above, but for all dots together

fAll = dMRCA_regression(L.Age,metric,L.Cutotype,clrs_all,MolecularClockRate,'all');

%% Stacked bars of all dMRCA

fStacked = dMRCA_stacked(L.Age,metric,L.Cutotype,clrs_all,MolecularClockRate,['all ' strrep(metricname,"_"," ") 'lineage per subject']);


end

%% runs plot subfunctions 
function p= dMRCA_regression(X,Yi,Cutotypes,clrs,MolecularClockRate,str)


p = figure;
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





%%