function f = TemporalDynamicsIndividual(A,ConcatenatedReadsDataSubjectTime,MinAssigned,cutoff_real,ttl)
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

mintps=3;

[uSID,~,iu] = unique(SID);
U=numel(uSID);
[AbundancesTPs,SidsTPs,TPsTPs,CutoTPs]=deal([]);
for u = 1:U
    MindexSubject = iu==u;
    UniqueTpsThisSid = unique(TP(MindexSubject));
    if numel(UniqueTpsThisSid)>=mintps
        N=numel(UniqueTpsThisSid);
        TPsTPs =[TPsTPs {(TP(MindexSubject))}];
        AbundancesTPs =[AbundancesTPs {A(MindexSubject,:)}];
        SidsTPs =[SidsTPs {repmat(unique(SID(MindexSubject)),1,N)}];
        CutoTPs=[CutoTPs {repmat(unique(Cutotype(MindexSubject)),1,N)}];
    end
end

%% differences between numbered timepoints

Tdiff=get_timepoint_differences;

%%
colors_all=[0.4000    0.7608    0.6471;0.9882    0.5529    0.3843; 0.5529    0.6275    0.7961; 0.9059    0.5412    0.7647; 0.6510    0.8471    0.3294];
CutoAll=string(cellfun(@unique,CutoTPs));
CutoAll(CutoAll=="3")="P";
uSID=cellfun(@unique,SidsTPs);
U=numel(uSID);

f=figure;
for i = 1:U
    Ntps = numel(TPsTPs{i});
    Nint=Ntps-1;
    IntervalDuration=arrayfun(@(x) {Tdiff(TPsTPs{i}(x),TPsTPs{i}(x+1))}, 1:(Nint));
    IntervalDuration=cellfun(@(x) {round(x,2,'decimals')},IntervalDuration);


    % Xdata is time from initial
    XTime = cumsum([0 horzcat(IntervalDuration{:})]);
    % Resize abundance only for lineages ever found
    A = AbundancesTPs{i};
    EverFound=sum(A>0,1)>0;
    A=A(:,EverFound);
    % Get colors, one for each clade
    N = size(A,2);
    clrs = cbrewer2('Paired',N);
    subplot(2,U,i)
    hold on
    title(([uSID(i)])) % ttl]))
    brs=bar(A,'stacked','FaceColor','flat');
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
    N = size(A,2);
    arrayfun(@(x) plot(XTime,A(:,x),'Color',clrs(x,:),'LineWidth',2), 1:N);
    xticks(XTime)
    ylim([0 1])
    xlabel('Time since T1')
    ylabel('Abundance')

end
