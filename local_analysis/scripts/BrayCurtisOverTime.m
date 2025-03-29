function T = BrayCurtisOverTime(A,M,MinAssigned,MinReal,ttl,min_interval_years,ManuscriptColors)

%% filter abundances
enough = sum(A,2)>=MinAssigned;

A=A(enough,:);
M=M(enough,:);

A(A<MinReal)=0;


%% parse needed variables from M


TP=M.TP;
SID=M.SID;
Cutotype=M.subject_cutotype;
Cutotype(contains(SID,"P"))=3;
Age=M.Age;

%% get differences between timepoints in years
Tdiff=get_timepoint_differences;

%% loop through subjects
[uSID,~,iu] = unique(SID);
U=numel(uSID);
[AbundancesTPs,SidsTPs,TPsTPs,CutoTPs,AgesTPs]=deal([]);
for u = 1:U
    MindexSubject = iu==u;
    UniqueTpsThisSid = unique(TP(MindexSubject));
    TPsubject = TP(MindexSubject);
    TimeSinceInitial = Tdiff(TPsubject(1),TPsubject(:));
    if numel(UniqueTpsThisSid)>=2&max(TimeSinceInitial)>min_interval_years
        N=numel(UniqueTpsThisSid);
        TPsTPs =[TPsTPs {(TP(MindexSubject))}];
        AbundancesTPs =[AbundancesTPs {A(MindexSubject,:)}];
        SidsTPs =[SidsTPs {repmat(unique(SID(MindexSubject)),1,N)}];
        CutoTPs=[CutoTPs {repmat(unique(Cutotype(MindexSubject)),1,N)}];
        AgesTPs=[AgesTPs {repmat(unique(Age(MindexSubject)),1,N)}];
        
    end
end


%% get color maps and measure dissimiliarity from self over time

clrs_all = ManuscriptColors.CutotypeColors;
TimeSinceInitial = cellfun(@(x) {Tdiff(x(1),x(:))},  TPsTPs);
BetaDissimilarityFromInitial = cellfun(@(x) {bray_curtis(x)},  AbundancesTPs);
BetaDissimilarityFromInitial=cellfun(@(x) {x(1,1:end)}, BetaDissimilarityFromInitial);


%% make plot


h=gca;
plot([TimeSinceInitial{1}],[BetaDissimilarityFromInitial{1}],'LineWidth',2,'Color',clrs_all(CutoTPs{1}(1),:))
hold on
title(ttl)
%arrayfun(@(x) plot([TimeSinceInitial{x}],[BetaDissimilarityFromInitial{x}],'LineWidth',1,'Color',clrs_all(CutoTPs{x}(1),:)), 1:numel(TimeSinceInitial))
arrayfun(@(x) plot([TimeSinceInitial{x}],[BetaDissimilarityFromInitial{x}],'LineWidth',2,'Color',clrs_all(CutoTPs{x}(1),:)), 2:numel(TimeSinceInitial))

yticks([0:.25:1])
ylim([0 1])
% yticks([-1:.5:2])
xlim([0 1.5])
xticks([0:.25:1.5])
% arrayfun(@(x) textscatter([max(TimeSinceInitial{x})],[BetaDissimilarityFromInitial{x}(end)],SidsTPs{x}(1)), 1:numel(TimeSinceInitial))
xlabel('Time')
ylabel('B-C Dissimilarity from T1')


%% fill table

T=table;
T.TimeSinceInitial=TimeSinceInitial';
T.BetaDissimilarityFromInitial=BetaDissimilarityFromInitial';
T.TPsTPs=TPsTPs';
T.SidsTPs=SidsTPs';
T.CutoTPs=CutoTPs';
T.Ys_simmilarity=BetaDissimilarityFromInitial';
T.AgeTPs=AgesTPs';
T.Abundances=AbundancesTPs';
%%