function clrs = getGroupProgressiveColors(Ages,ManuscriptColors)

Groups = ages2groups(Ages);
Groups(Groups==2)=1;
Groups(Groups==3)=2;

N = numel(Ages);
clrs = repmat(ManuscriptColors.CutotypeColors(3,:),N,1);


AgeRange = floor(min(Ages(Groups==1))):ceil(max(Ages(Groups==1)));

Ucolors = cbrewer2(ManuscriptColors.AgeColorMap,numel(AgeRange));


Group1Colors = arrayfun(@(x) {Ucolors(AgeRange==x,:)}, round(Ages(Groups==1)));

clrs(Groups==1,:)=vertcat(Group1Colors{:});

figure;
scatter(Ages,ones(size(Ages)),'filled','CData',clrs)