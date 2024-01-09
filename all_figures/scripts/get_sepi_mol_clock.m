%aggregate all

function [rate,fg]= get_sepi_mol_clock(SubjectClusterPairs,tps,root_to_tip)
%% for all pairs, compare root-to-tip distance over time 

T1="2018-05-29";
T2="2018-10-25";
T3="2018-12-13";
T4="2019-06-04";
T5="2019-10-24";
T6="2022-12-02";

Tall = [T1 T2 T3 T4 T5 T6];
Tall=datetime(Tall);
Tdiff = abs(years(Tall-Tall'));
%

Yall = [];
Xall=[];
is_sepi=[];
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
        is_sepi=[is_sepi repmat(contains(SubjectClusterPairs(i),"sepi"),1,numel(ydeltas))];
    end
end

%
X=Xall(is_sepi==1);
Y=Yall(is_sepi==1);

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
text(0,13,['y=' char(string(rate)) 'x']);
text(0,12,['adjR-squared=' char(string(adjr2))]);
text(0,11,['pvalue=' char(string(Pval))]);