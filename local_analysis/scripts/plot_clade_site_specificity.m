function [f,pvs_all] = plot_clade_site_specificity(M,A,str)
%% initialize figure
f=figure;


%% remove comparisons where <70% of community is assigned
Assigned = sum(A,2);
Enough = Assigned>.7;
A(A<.01)=0;


%% get bray-curtis dissimilarity 
BC = bray_curtis(A);
BC(~(Enough&Enough'))=nan;% remove instances where both samples being compared don't have enough assigned


%% parse which sites, subject, timepoints, and families corrospond to each
[sites1, sites2] = deal(M.location);sites2=sites2';
[SID1,SID2]=deal(M.SID);SID2=SID2';
[TP1, TP2] = deal(M.TP);TP2=TP2';
[Fam1,Fam2]=deal(get_family_numbers(M.SID));Fam2=Fam2';

% make all possible site-specific comparison booleans
FF = sites1=="F"&sites2=="F";
NN = sites1=="N"&sites2=="N";
KK = sites1=="K"&sites2=="K";
CC = sites1=="C"&sites2=="C";
FN = (sites1=="F"&sites2=="N")|(sites1=="F"&sites2=="N")';
FK = (sites1=="F"&sites2=="K")|(sites1=="F"&sites2=="K")';
FC = (sites1=="F"&sites2=="C")|(sites1=="F"&sites2=="C")';
NK = (sites1=="N"&sites2=="K")|(sites1=="N"&sites2=="K")';
NC = (sites1=="N"&sites2=="C")|(sites1=="N"&sites2=="C")';
KC = (sites1=="K"&sites2=="C")|(sites1=="K"&sites2=="C")';

FOther=(sites1=="F"&sites2~="F")|(sites2=="F"&sites1~="F");
Nother=(sites1=="N"&sites2~="N")|(sites2=="N"&sites1~="N");
Kother=(sites1=="K"&sites2~="K")|(sites2=="K"&sites1~="K");
Cother=(sites1=="C"&sites2~="C")|(sites2=="C"&sites1~="C");

% booleans which are true when comparing same subject, timepoint, famility, etc
SameTimepoint = TP1==TP2;
SameSubject = SID1==SID2;
SameFamily=Fam1==Fam2;
upoint = tril(ones(size(SameTimepoint)),-1)==1;% the lower triangle of the matrix, such that you aren't including duplicate points


%% Pairwise comparison of sites between unrelated people
Ys = [{BC(~SameFamily&FF&upoint)} {BC(upoint&~SameFamily&FOther)} {BC(~SameFamily&NN&upoint)} {BC(upoint&~SameFamily&Nother)} {BC(~SameFamily&KK&upoint)} {BC(upoint&~SameFamily&Kother)} {BC(upoint&~SameFamily&CC&upoint)} {BC(~SameFamily&Cother)}];
Xs = arrayfun(@(x) {x*ones(size(Ys{x}))}, 1:numel(Ys));

% concatenate together for plot
Ys= vertcat(Ys{:});
Xs = vertcat(Xs{:});

R = arrayfun(@(x) {randsample(find(Xs==x),1000)}, unique(Xs));
R =vertcat(R{:});
% takes forever to plot, so take a subsample
XP = Xs(R);
YP = Ys(R);

% make plot, all gray
subplot(2,1,1)

swarmchart(XP,YP,5,'filled','MarkerFaceColor',[.85 .85 .85])
hold on

% add median lines and p-values
for x = 1:8
    mn = median(Ys(x==Xs),'omitnan');
    plot([x-.25 x+.25],[mn mn],'Color','black','LineWidth',1)
end

% get p-values for comparisons
ps= [0 0 0 0];
idx=1;
for i=[1 3 5 7]
     [ps(idx)]=ranksum(Ys(Xs==i),Ys(Xs==(1+i)));

    idx=idx+1;
end
ps = round(ps,2,'significant');
ps = string(ps);

% add p-values to plot
idx=1:2:8;
arrayfun(@(x) text(idx(x)+.3,1.06,['p=' char(ps(x))]), 1:4)
arrayfun(@(x) plot([x x+1],[1.05 1.05],'Color','black'), 1:2:8)

% add title, ticks, and labels
title([ str  ' Pairwise dissimilarity across unrelated people'])
xticks([1:8])
xticklabels(["Fh-Fh" "Fh-Other" "Ns-Ns" "Ns-Other" "Ck-Ck" "Ck-Other" "Ch-Ch" "Ch-Other"])
ylabel("B-C Dissimilarity")

yticks([0:.2:1])

%% Compare across-site differences on the same subject vs unrelated people
Ys = [{BC(SameSubject&SameTimepoint&FOther)} {BC(~SameFamily&FOther)} ... 
    {BC(SameSubject&SameTimepoint&Nother)} {BC(~SameFamily&Nother)} ... 
    {BC(SameSubject&SameTimepoint&Kother)} {BC(~SameFamily&Kother)}... 
    {BC(SameSubject&SameTimepoint&Cother)} {BC(~SameFamily&Cother)}];

Xs = arrayfun(@(x) {x*ones(size(Ys{x}))}, 1:numel(Ys));

% concatenate for plot
Xs=vertcat(Xs{:});
Ys=vertcat(Ys{:});

clrs = repmat([.85 .85 .85],numel(Ys),1);
odd = mod(Xs,2)==1;
clrs(odd,:)=repmat([0 1 1],sum(odd),1);

R = arrayfun(@(x) {randsample(find(Xs==x),min([1000 sum(Xs==x)]))}, unique(Xs));
R =vertcat(R{:});
% takes forever to plot, so take a subsample
XP = Xs(R);
YP = Ys(R);
clrs=clrs(R,:);

% make plot, all gray
subplot(2,1,2)
swarmchart(XP,YP,5,'filled','CData',clrs)
hold on

% add median lines and p-values
for x = 1:8
    mn = median(Ys(x==Xs),'omitnan');
    plot([x-.25 x+.25],[mn mn],'Color','black','LineWidth',1)
end

% get p-values for comparisons
ps2= [0 0 0 0];
idx=1;
for i=[1 3 5 7]
          [ps2(idx)]=ranksum(Ys(Xs==i),Ys(Xs==(1+i)));

    idx=idx+1;
end
ps2 = round(ps2,2,'significant');
ps2 = string(ps2);

% add p-values to plot
idx=1:2:8;
arrayfun(@(x) text(idx(x)+.3,1.06,['p=' char(ps2(x))]), 1:4)
arrayfun(@(x) plot([x x+1],[1.05 1.05],'Color','black'), 1:2:8)

xticks([1.5:2:9])
xticklabels(["Fh-Other" "Ns-Other" "Ck-Other" "Ch-Other"])

text(8.5,1,"SameSubject,Same Timepoint",'Color',[0 1 1]);
text(8.5,.85,"UnrelatedSubject",'Color',[.85 .85 .85]);
ylabel("B-C Dissimilarity")
yticks([0:.2:1])
pvs_all = [ps ps2]