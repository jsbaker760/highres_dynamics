function [f1,f2,f3]=plot_molecularclocks_individual(SubjectClusterPairs,tps,root_to_tip,MsubjectTime)

%% number of instances to iterate
C=numel(SubjectClusterPairs);

%% Initialize correlation coefficients and

[rho,pval]=deal(zeros(numel(root_to_tip),1));
hascore = false(size(rho));
[TimeRottToTipAll,AllRootToTipAllAll] = deal(cell(C,1));

% iterate
for c = 1:C
    if ~isempty(root_to_tip{c})
        SubjectClusterPairs(c)
        timepoints = tps{c};
        Xs = arrayfun(@(x) x*ones(numel(root_to_tip{c}{timepoints==x}),1) , timepoints , 'UniformOutput' , false);
        Xs = vertcat(Xs{:});
        Ys =horzcat(root_to_tip{c}{:});
        [rho(c),pval(c)] = corr(Xs,Ys');
        TimeRottToTipAll{c}=Xs;
        AllRootToTipAllAll{c}=Ys;
        hascore(c)=true;
    end
end

%% booleans for both species to split data

CacnesBool = contains(SubjectClusterPairs,"cacnes")&hascore';
SepiBool = contains(SubjectClusterPairs,"sepi")&hascore';

%% C. acnes corr coeficient vs p-value with FDR cutoff

f1=figure;hold on;
scatter(rho(CacnesBool),log10(pval(CacnesBool)),'filled','MarkerFaceColor',[.85 .85 .85])
xlim([-1 1])

FDR = .01;
pvs = pval(CacnesBool);
crit_val = ((1:numel(pvs))./numel(pvs)).*FDR;

[ps_sorted]= sort(pvs(:));

crit_p=max(ps_sorted(ps_sorted<crit_val'));

significant_p = pvs<crit_p;

plot([-1 1],[log10(crit_p) log10(crit_p)],'LineStyle',':','Color','red','LineWidth',2)

rhosig = rho(CacnesBool);rhosig=rhosig(significant_p);
pvalsig = pval(CacnesBool);pvalsig=pvalsig(significant_p);
scatter(rhosig,log10(pvalsig),'filled','MarkerFaceColor','green')
H=gca;
H.LineWidth=2;
xticks([-1:.5:1])
ylim([-8 0])
text(-1,log10(crit_p)+.2,['crit=' char(string(crit_p))])
title(['C. acnes clock values, FDR=' char(string(FDR))])
xlabel('Pearson correlation coefficient')
ylabel('log10 p-value')

%% As above for S. epidermidis

f2=figure;hold on;
scatter(rho(SepiBool),log10(pval(SepiBool)),'filled','MarkerFaceColor',[.85 .85 .85])
xlim([-1 1])

FDR = .01;
pvs = pval(SepiBool);
crit_val = ((1:numel(pvs))./numel(pvs)).*FDR;

[ps_sorted]= sort(pvs(:));

crit_p=max(ps_sorted(ps_sorted<crit_val'));

significant_p = pvs<=crit_p;

plot([-1 1],[log10(crit_p) log10(crit_p)],'LineStyle',':','Color','red','LineWidth',2)

rhosig = rho(SepiBool);rhosig=rhosig(significant_p);
pvalsig = pval(SepiBool);pvalsig=pvalsig(significant_p);
scatter(rhosig,log10(pvalsig),'filled','MarkerFaceColor','green')
H=gca;
H.LineWidth=2;
xticks([-1:.5:1])
ylim([-8 0])
text(-1,log10(crit_p)+.2,['crit=' char(string(crit_p))])
title(['S. epi. clock values, FDR=' char(string(FDR))])
xlabel('Pearson correlation coefficient')
ylabel('log10 p-value')

%%

SignificantSepiInstances = SubjectClusterPairs(SepiBool);
SignificantSepiInstances=SignificantSepiInstances(significant_p);

%%

f3=figure;
SepiCoeffs = zeros(numel(SignificantSepiInstances),2);
N=numel(SignificantSepiInstances);
for i = 1:N
    subplot(N,1,i)
    Xs = TimeRottToTipAll{SubjectClusterPairs==SignificantSepiInstances(i)};
    Ys = AllRootToTipAllAll{SubjectClusterPairs==SignificantSepiInstances(i)};
    Ym=ceil(max(Ys));
    Ymn=floor(min(Ys))
    Subj=strsplit(SignificantSepiInstances(i)," ");
    Ages = arrayfun(@(x) MsubjectTime.Age(find(MsubjectTime.SID==Subj(1)&MsubjectTime.TP==x,1)), Xs);
    Xs=Ages;
    swarmchart(Xs,Ys,'filled','MarkerFaceColor','black','XJitter','density','XJitterWidth',.2)
    hold on
    SepiCoeffs(i,:)=coeffvalues(fit(Xs,Ys','poly1'));
    coeff=corr(Xs,Ys',"type","Spearman")
    Xs = 1:Ages+10;
    Ys = (SepiCoeffs(i,1)*Xs)+SepiCoeffs(i,2);
    plot(Xs,Ys,'LineWidth',1,'Color','black','LineStyle','--')
    h=gca;
    h.LineWidth=2;  
    title(SignificantSepiInstances{i})
    xlim([floor(min(Ages)) ceil(max(Ages))])
    ylim([Ymn Ym])
    eq= ['f(x)=' char(string(round(SepiCoeffs(i,1),2,'significant'))) 'x+' char(string(round(SepiCoeffs(i,2),2,'significant')))];
    R2 = char(string(round(coeff^2,2,"significant")))
    text(floor(min(Ages)),Ym-1,eq)
    text(floor(min(Ages)),Ym-3,['R^2=' R2])
    pbaspect([1 1 1])

end
% 
