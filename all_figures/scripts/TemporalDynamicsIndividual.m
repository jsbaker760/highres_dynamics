function T = TemporalDynamicsIndividual(A,ConcatenatedReadsDataSubjectTime,MinAssigned,cutoff_real,ttl)
%% First, normalize abundances and remove lowly abundant clades

enough = sum(A,2)>=MinAssigned;
A=A(enough,:);
TP=ConcatenatedReadsDataSubjectTime.TP(enough);
SID=ConcatenatedReadsDataSubjectTime.SID(enough);
Cutotype=ConcatenatedReadsDataSubjectTime.subject_cutotype(enough);
Cutotype(contains(SID,"P"))=3;

A(sum(A,2)>1,:)=A(sum(A,2)>1,:)./sum(A(sum(A,2)>1,:),2);
A(A<cutoff_real)=0;
Age=ConcatenatedReadsDataSubjectTime.Age(enough);
%%
load data/ColorBlind.mat
%%
mintps=3;

%%

[uSID,~,iu] = unique(SID);
U=numel(uSID);
[AbundancesTPs,SidsTPs,TPsTPs,CutoTPs,AgesTPs]=deal([]);
for u = 1:U
    MindexSubject = iu==u;
    UniqueTpsThisSid = unique(TP(MindexSubject));
    if numel(UniqueTpsThisSid)>=mintps
        N=numel(UniqueTpsThisSid);
        TPsTPs =[TPsTPs {(TP(MindexSubject))}];
        AbundancesTPs =[AbundancesTPs {A(MindexSubject,:)}];
        SidsTPs =[SidsTPs {repmat(unique(SID(MindexSubject)),1,N)}];
        CutoTPs=[CutoTPs {repmat(unique(Cutotype(MindexSubject)),1,N)}];
        AgesTPs=[AgesTPs {repmat(unique(Age(MindexSubject)),1,N)}];
        if numel(unique(Age(MindexSubject)))>1
            foo=1
        end
    end
end

%%


Tdiff=get_timepoint_differences;

%%

clrs_all = [ColorBlind.LightBlue; ColorBlind.Blue; ColorBlind.Vermillion];
TimeSinceInitial = cellfun(@(x) {Tdiff(x(:),x(:))},  TPsTPs);
BetaDissimilarityFromInitial = cellfun(@(x) {jaccard(x)},  AbundancesTPs);


%%
colors_all=[0.4000    0.7608    0.6471;0.9882    0.5529    0.3843; 0.5529    0.6275    0.7961; 0.9059    0.5412    0.7647; 0.6510    0.8471    0.3294];
CutoAll=string(cellfun(@unique,CutoTPs))
CutoAll(CutoAll=="3")="P"
% F=factor(U);
% I=min(F);
% J=U/I;
uSID=cellfun(@unique,SidsTPs)
U=numel(uSID)

% F = figure;

[IntervalSIDs,IntervalGains,IntervalLoss,IntervalStart,IntervalEnd,IntervalAge,IntervalCutotypes,IntervalTime,IntervalBeta]=deal([]);
for i = 1:U
     
    Ntps = numel(TPsTPs{i});
    Nint=Ntps-1;
    IntervalXY=arrayfun(@(x) {AbundancesTPs{i}([x x+1],:)}, 1:(Nint));
    IntervalDuration=arrayfun(@(x) {Tdiff(TPsTPs{i}(x),TPsTPs{i}(x+1))}, 1:(Nint));
    IntervalDuration=cellfun(@(x) {round(x,2,'decimals')},IntervalDuration);
    IntervalBC=arrayfun(@(x) {BetaDissimilarityFromInitial{i}(x,x+1)}, 1:(Nint));
    AbundancesGained=cellfun(@(x) {x(2,x(1,:)==0&x(2,:)>0)}, IntervalXY);
    AbundancesLost=cellfun(@(x) {x(1,x(1,:)>00&x(2,:)==0)}, IntervalXY);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IntervalSIDs=[IntervalSIDs repmat(uSID(i),1,Nint)];
    IntervalGains=[IntervalGains AbundancesGained];
    IntervalLoss=[IntervalLoss AbundancesLost];
    IntervalStart=[IntervalStart     cellfun(@(x) {x(1,:)}, IntervalXY)];
    IntervalEnd=[IntervalEnd     cellfun(@(x) {x(2,:)}, IntervalXY)];
    IntervalAge=[IntervalAge {repmat(unique(AgesTPs{i}),1,Nint)}];
    IntervalCutotypes=[IntervalCutotypes repmat(unique(CutoTPs{i}),1,Nint)];
    IntervalTime=[IntervalTime IntervalDuration];
    IntervalBeta=[IntervalBeta IntervalBC];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(U,3,((i-1)*3)+1)
%     subplot(1,3,1)
%     hold on
%     arrayfun(@(x) scatter(IntervalXY{x}(1,:),IntervalXY{x}(2,:),'filled','MarkerFaceColor',colors_all(x,:)), 1:Nint)
%     xlabel('Abundance T1')
%     ylabel('Abundance T2')
%     xlim([0 1])
%     ylim([0 1])
% %     l=legend;
% %     l.String=arrayfun(@(x) strjoin(["Interval" string(x) "Duration" string(round(IntervalDuration{x},2)) "years"]), 1:Nint);
%     pbaspect([1 1 1])
% 
% %     subplot(U,3,((i-1)*3)+2)
%     subplot(1,3,2)
%     XsGained=1:2:(((Nint-1)*2)+1)
%     XsLost=2:2:(((Nint-1)*2)+2)
%     Xs=[XsGained XsLost];
%     Ys=[AbundancesGained AbundancesLost]
%     lbls= repmat(["Gained" "Lost"],1,Nint)
%     Xs=arrayfun(@(x) {repmat(Xs(x),1,numel(Ys{x}))}, 1:numel(Xs));
%     Xs=horzcat(Xs{:})
%     colors_int = arrayfun(@(x) {repmat(colors_all(x,:),2,1)}, 1:Nint)
%     colors_int=vertcat(colors_int{:});
%     colors_int=colors_int(Xs,:);
%     Ys=horzcat(Ys{:})
%     swarmchart(Xs,Ys,'filled','cdata',colors_int)
%     xticks(1:max(XsLost))
%     xticklabels(lbls)
%     xlim([0 max(XsLost+1)])
%     ylim([0 1])
%     pbaspect([1 1 1])
%     ylabel('Abundance')
% 
% %     subplot(U,3,((i-1)*3)+3)
%     subplot(1,3,3)
%     hold on
%     scatter(horzcat(IntervalDuration{:}),horzcat(IntervalBC{:}),'filled','CData',colors_all(1:Nint,:))
%     xlim([0 1]);
%     ylim([0 1])
%     xlabel('Interval Duration')
%     ylabel('B-C Dissimilarity')
%     pbaspect([1 1 1])
%     sgtitle(strjoin([ "Subject" uSID(i) "Cutotype" CutoAll(i) ttl]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Xdata is time from initial
    XTime = cumsum([0 horzcat(IntervalDuration{:})]);
    % Resize abundance only for lineages ever found
    A = AbundancesTPs{i};
    EverFound=sum(A>0,1)>0;
    A=A(:,EverFound);
    % Get colors, one for each clade
    N = size(A,2);
    clrs = cbrewer2('Paired',N);
    % Plot data as stacked bars
%     f=figure('Renderer', 'painters', 'Position', [680,616,231,361])
    subplot(2,U,i)
    hold on
    title(([uSID(i)])) % ttl]))
    brs=bar(A,'stacked','FaceColor','flat')
    for b = 1:numel(brs)
        brs(b).CData=clrs(b,:);
    end
    ylim([0 1])
    xticks([1:size(A,1)])
    xticklabels(TPsTPs{i})
    % plot abundance of inidividual lineages
    xlabel('Timepoints')
    ylabel('Abundance')
    subplot(2,U,i+U)
    hold on
%     Stable=sum(A>0,1)==size(A,1);
%     A = A(:,~Stable);
%     clrs = clrs(~Stable,:);
    N = size(A,2);
    arrayfun(@(x) plot(XTime,A(:,x),'Color',clrs(x,:),'LineWidth',2), 1:N);
    xticks(XTime)
    
    ylim([0 1])
    xlabel('Time since T1')
    ylabel('Abundance')
%     saveas(f,['ManuscriptFigures/' char(strjoin([uSID(i) ttl])) '.svg'],'svg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

T=table;
T.IntervalSIDs=IntervalSIDs';
T.IntervalGains=IntervalGains';
T.IntervalLoss=IntervalLoss';
T.IntervalStart=IntervalStart';
T.IntervalEnd=IntervalEnd';
% T.IntervalAge=IntervalAge';
T.IntervalCutotypes=IntervalCutotypes';
T.IntervalTime=horzcat(IntervalTime{:})';
T.IntervalBeta=horzcat(IntervalBeta{:})';
