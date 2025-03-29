function [fMin,fMax,fAll]=compare_lineage_dmrcas_all(L,MolecularClockRate,metric,clrsgroups)

%% average distances
metric = cellfun(@mean,metric);

%% get unique subjects and their cutotypes

[uSID,uidx] = unique(L.SID);
Cutotypes = L.Cutotype(uidx);

%% colors for each dot
clrs_subjects=clrsgroups(Cutotypes,:);


%% get the age of subject (X) and the mean, min, max Y values per subject

X = arrayfun(@(x) mean(L.Age(x==L.SID)), uSID);
Ymin = arrayfun(@(x) min(metric(x==L.SID)), uSID);
Ymax = arrayfun(@(x) max(metric(x==L.SID)), uSID);

%% plot min, max, and all dMRCA with one dot per subject

    clrs_all=clrsgroups(L.Cutotype,:);
    fAll = dMRCA_regression(L.Age,metric,L.Cutotype,clrs_all,MolecularClockRate,'all');
    fMin=dMRCA_regression(X,Ymin,Cutotypes,clrs_subjects,MolecularClockRate,['min']);
    fMax=dMRCA_regression(X,Ymax,Cutotypes,clrs_subjects,MolecularClockRate,['max']);

end