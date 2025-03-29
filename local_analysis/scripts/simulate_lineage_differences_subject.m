function [f, SimData]=simulate_lineage_differences_subject(SID,A,DM,T,replacement)
%% remove rows without enough data
enough=sum(A,2,'omitnan')>0;
A=A(enough,:);

SID=SID(enough);

[~,~,ic]=unique(SID);
%% get dimensions and incides for permutations
[S,C]=size(A);

%% get indices to permute depending on if you decide with or without replacement
if replacement
    perms = randi(C, T, C);
else
    [~,perms] = sort(rand(T,C),2);
end

%% since unweighted, always equals 1 (boolean is also faster)
A=A>0;

%% set middle of DM to NaN to prevent calculation of self-to-self
DM(eye(size(DM))==1)=NaN;

%% Reshape permutations
permuted_A=reshape(A(:,perms'),S,C,T);

%% get real values for mean, min, max pairwise distance to compare to simulationss
AverageWithinClusteDistances=NaN(S,1);
MinWithinClusteDistances=NaN(S,1);
MaxWithinClusteDistances=NaN(S,1);
for s=1:S
    b = A(s,:);
    if sum(b)>1
        DMs=DM(b,b);
        AverageWithinClusteDistances(s)=mean(DMs(:),'omitnan');
        MinWithinClusteDistances(s)=min(DMs(:),[],1,'omitnan');
        MaxWithinClusteDistances(s)=max(DMs(:),[],1,'omitnan');
    end
end


AverageWithinClusteDistances=arrayfun(@(x) nanmean(AverageWithinClusteDistances(ic==x)), 1:max(ic));
MinWithinClusteDistances=arrayfun(@(x) nanmean(MinWithinClusteDistances(ic==x)), 1:max(ic));
MaxWithinClusteDistances=arrayfun(@(x) nanmean(MaxWithinClusteDistances(ic==x)), 1:max(ic));

%% Run simulations
SimulatedAverageWithinClusteDistances=NaN(T,S);
SimulatedMinWithinClusteDistances=NaN(T,S);
SimulatedMaxWithinClusteDistances=NaN(T,S);

for t=1:T
    % simulation instance
    for s=1:S
        b=permuted_A(s,:,t);
        if sum(b)>1
            DMs=DM(b,b);
            SimulatedAverageWithinClusteDistances(t,s)=mean(DMs(:),'omitnan');
            SimulatedMinWithinClusteDistances(t,s)=min(DMs(:),[],1,'omitnan');
            SimulatedMaxWithinClusteDistances(t,s)=max(DMs(:),[],1,'omitnan');

        end
    end
end

SimulatedAverageWithinClusteDistances=arrayfun(@(x) {nanmean(SimulatedAverageWithinClusteDistances(:,ic==x),2)}, 1:max(ic));
SimulatedMinWithinClusteDistances=arrayfun(@(x) {nanmean(SimulatedMinWithinClusteDistances(:,ic==x),2)}, 1:max(ic));
SimulatedMaxWithinClusteDistances=arrayfun(@(x) {nanmean(SimulatedMaxWithinClusteDistances(:,ic==x),2)}, 1:max(ic));

SimulatedAverageWithinClusteDistances=horzcat(SimulatedAverageWithinClusteDistances{:});
SimulatedMinWithinClusteDistances=horzcat(SimulatedMinWithinClusteDistances{:});
SimulatedMaxWithinClusteDistances=horzcat(SimulatedMaxWithinClusteDistances{:});

%% avergae across simulation
SimulatedAverageWithinClusteDistances=mean(SimulatedAverageWithinClusteDistances,2,'omitnan');
SimulatedMinWithinClusteDistances=mean(SimulatedMinWithinClusteDistances,2,'omitnan');
SimulatedMaxWithinClusteDistances=mean(SimulatedMaxWithinClusteDistances,2,'omitnan');

%% average of real data
MeanReal=mean(AverageWithinClusteDistances,'omitnan');
MinReal=mean(MinWithinClusteDistances,'omitnan');
MaxReal=mean(MaxWithinClusteDistances,'omitnan');



%% Create structure to add simulation data
SimData=struct;

%% make a figure
f=figure;

% Average within cluster distances
bw=100;
Dsim=SimulatedAverageWithinClusteDistances;
Dreal=MeanReal;
subplot(1,3,1)
h=histogram(Dsim,1000,'FaceColor',[.85 .85 .85],'EdgeAlpha',0,'Normalization','probability');
hold on
ylim([0 4e-3])
yticks([0:1e-3:4e-3])
plot([Dreal Dreal],[0 4e-3],'Color','red','LineWidth',2)
% obs is either the observed value or 1-obs to get 2-tail result
obs=Dreal;
pv=mean(Dsim>=obs);
pv=min([pv 1-pv]);
pv = pv*2;
SimData.p_mean=pv;
pv=['p=' char(string(round(pv,2,'significant')))];
text(h.BinLimits(1),max(h.Values),pv)
% label
xlabel('average inter-cluster distance (averaged)')
ylabel('proportion simulations')
pbaspect([1 1 1])

% Min within cluster distances
Dsim=SimulatedMinWithinClusteDistances;
Dreal=MinReal;
subplot(1,3,2)
h=histogram(Dsim,1000,'FaceColor',[.85 .85 .85],'EdgeAlpha',0,'Normalization','probability');
hold on
ylim([0 4e-3])
yticks([0:1e-3:4e-3])
plot([Dreal Dreal],[0 4e-3],'Color','red','LineWidth',2)
% obs is either the observed value or 1-obs to get 2-tail result
obs=Dreal;
pv=mean(Dsim>=obs);
pv=min([pv 1-pv]);
pv = pv*2;
SimData.p_min=pv;
pv=['p=' char(string(round(pv,2,'significant')))];
text(h.BinLimits(1),max(h.Values),pv)
% label
xlabel('minimum inter-cluster distance (averaged)')
ylabel('proportion simulations')
pbaspect([1 1 1])

% Max within cluster distances
Dsim=SimulatedMaxWithinClusteDistances;
Dreal=MaxReal;
subplot(1,3,3)
h=histogram(Dsim,1000,'FaceColor',[.85 .85 .85],'EdgeAlpha',0,'Normalization','probability');
hold on
ylim([0 4e-3])
yticks([0:1e-3:4e-3])
plot([Dreal Dreal],[0 4e-3],'Color','red','LineWidth',2)% obs is either the observed value or 1-obs to get 2-tail result
obs=Dreal;
pv=mean(Dsim>=obs);
pv=min([pv 1-pv]);
pv = pv*2;
SimData.p_max=pv;
pv=['p=' char(string(round(pv,2,'significant')))];
text(h.BinLimits(1),max(h.Values),pv)
% label
xlabel('maximum inter-cluster distance (averaged)')
ylabel('proportion simulations')
pbaspect([1 1 1])

%% Add simulations to data structure
SimData.SimMax=SimulatedMaxWithinClusteDistances;
SimData.SimMean=SimulatedAverageWithinClusteDistances;
SimData.SimMin=SimulatedMinWithinClusteDistances;
SimData.MaxReal=MaxReal;
SimData.MeanReal=MeanReal;
SimData.MinReal=MinReal;
