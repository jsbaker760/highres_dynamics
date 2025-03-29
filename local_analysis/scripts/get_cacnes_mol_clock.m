function [rate,fg]= get_cacnes_mol_clock(SubjectClusterPairs,tps,root_to_tip)
%% for all pairs, compare root-to-tip distance over time 

Tdiff = get_timepoint_differences;
%

Yall = [];
Xall=[];
is_cacnes=[];
for i = 1:numel(SubjectClusterPairs)
    T = tps{i};
    if ~isempty(T)
        Neach = cellfun(@numel,root_to_tip{i});
        T1=T(1);
        xdeltas = arrayfun(@(x) Tdiff(T1, x), T);
        meanrtt=cellfun(@mean,root_to_tip{i});
        ydeltas = arrayfun(@(X) meanrtt(X)-meanrtt(X-1) , 2:(numel(meanrtt)));
        Yall = [Yall ydeltas];
        Xall = [Xall (xdeltas(2:end))'];
        is_cacnes=[is_cacnes repmat(contains(SubjectClusterPairs(i),"cacnes"),1,numel(ydeltas))];
    end
end

%
X=Xall(is_cacnes==1);
Y=Yall(is_cacnes==1);

%% initial figure

fg=figure;
scatter(X,Y,'filled')
hold on

%
ft1 = fittype({'x'});
[f,gof] = fit(X',Y',ft1);
rate=f.a;

plot(f,X,Y)

ylabel('mean(dMRCA)-D1')
xlabel('delta T')

%
[~,Pval]=corr(X',Y','type','Pearson');
adjr2=gof.adjrsquare;

xlim([0 1.1]);
text(0,3,['y=' char(string(rate)) 'x']);
text(0,2,['adjR-squared=' char(string(adjr2))]);
text(0,1,['pvalue=' char(string(Pval))]);