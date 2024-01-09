function plot_all_lineages_family(MsubjectTime)

%%

cacnes = MsubjectTime.CombinedCacnesLineages;
sepi=MsubjectTime.CombinedSepiLineages;

%%

Acacnes = sum(cacnes,2);
Asepi = sum(sepi,2);

MsubjectTime_cacnes = MsubjectTime(Acacnes>.7,:)
MsubjectTime_sepi = MsubjectTime(Asepi>.7,:)
cacnes=cacnes(Acacnes>.7,:);
sepi=sepi(Asepi>.7,:);
%%
family_1_sids = ["1AA" "1PA" "1PB"]

MsubjectTime_cacnes=MsubjectTime_cacnes(ismember(MsubjectTime_cacnes.SID,family_1_sids),:)
MsubjectTime_sepi=MsubjectTime_sepi(ismember(MsubjectTime_sepi.SID,family_1_sids),:)
cacnes=cacnes(ismember(MsubjectTime_cacnes.SID,family_1_sids),:)
sepi=sepi(ismember(MsubjectTime_sepi.SID,family_1_sids),:)
%
%%

all_cacnes_colors = sum(cacnes,1)>0
all_sepi_colors = sum(sepi,1)>0

sepi= sepi(:,all_sepi_colors);
cacnes= cacnes(:,all_cacnes_colors);

%%
cacnes_colors = cbrewer2('Spectral',sum(all_cacnes_colors))
sepi_colors = cbrewer2('PrGn',sum(all_sepi_colors))

%%

f = figure;
for i =1:3
subplot(3,1,i)
    bool=family_1_sids(i)==MsubjectTime_cacnes.SID;
    T = MsubjectTime_cacnes(bool,:);
    A = cacnes(bool,:);
    Xs = T.TP;

    title(family_1_sids(i))
    brs=bar(Xs,A,'stacked')
    for b = 1:16
        brs(b).FaceColor=cacnes_colors(b,:)
    end
    xticks([1:6])
    xlim([0 7])
    xlabel('sampling timepoint')
    ylabel('assigned abundance')
    title(family_1_sids(i))
end
sgtitle('C. acnes lineages in family 1')
%%

f = figure;
for i =1:3
subplot(3,1,i)
    bool=family_1_sids(i)==MsubjectTime_sepi.SID;
    T = MsubjectTime_sepi(bool,:);
    A = sepi(bool,:);
    Xs = T.TP;

    title(family_1_sids(i))
    brs=bar(Xs,A,'stacked')
    for b = 1:8
        brs(b).FaceColor=sepi_colors(b,:)
    end
    xticks([1:6])
    xlim([0 7])
    xlabel('sampling timepoint')
    ylabel('assigned abundance')
    title(family_1_sids(i))
end
sgtitle('S. epidermidis lineages in family 1')