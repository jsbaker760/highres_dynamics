function p=plot_cutotypes_per_tp(M,ManuscriptColors)

% get IDs of each subject (uSubject)

uSubect=unique(M.SID);
U=numel(uSubect);
uTPs = unique(M.TP);
T = numel(uTPs);

%% make X and Y data

% all possible TPs a subject could have been sampled At

% X is the number of the row to which the subject belongs in uSubject
Xs = repmat(1:numel(uSubect),T,1);
% Y is timepoint
Ys = repmat(uTPs',1,U);

%% Fill in data array to color dots by cutotype

D = numel(Ys);

% fills in UMAP cluster for each element for subject X and TP Y, 
% and is zero where the subject was not sampled at that TP
% then resizes it to be the same size as Y
SubjectTimeCutotype = reshape(arrayfun(@(x) max([M.UMAPClustersSubjectTime(M.SID==uSubect(Xs(x))&M.TP==Ys(x)) 0]), 1:D),T,U);

%% resize data to be 1D

Xs=Xs(:);
Ys=Ys(:);
C = SubjectTimeCutotype(:);

%% make colormap

clrs = repmat({[1 1 1]},D,1);
clrs(C==1)={ManuscriptColors.CutotypeColors(1,:)};
clrs(C==2)={ManuscriptColors.CutotypeColors(2,:)};
clrs = vertcat(clrs{:});

%% plot

p = figure;

scatter(Xs,Ys,100,'filled','CData',clrs,'Marker','square')

pbaspect([U T 1]);

xticks([1:U])
xticklabels(uSubect)
yticks(1:T)
yticklabels(1:T)
xlabel('Subject')
ylabel('Timepoint')



