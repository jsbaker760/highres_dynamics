function clrs = getGroupProgressiveColors(Ages,ManuscriptColors)

%%
%{

N = numel(Ages);

Groups = ages2groups(Ages);
Groups(Groups==2)=1;
Groups(Groups==3)=2;


MinAge= min(Ages(Groups==1))
MaxAge= max(Ages(Groups==1))

MaxBlue = .8;
MinBlue=0;

DiffColors = abs(MaxBlue-MinBlue)
DiffAge = abs(MaxAge-MinAge)

BlueAll = arrayfun(@(x) MinBlue + DiffColors*(abs(x-MinAge)/DiffAge), Ages);

%
clrs = zeros(N,3);
clrs(Groups==2,:)=.5;
clrs(Groups==1,3)=1
clrs(Groups==1,2)=BlueAll(Groups==1)
clrs(Groups==1,1)=BlueAll(Groups==1)

UcolorsChild = unique(clrs(Groups==1,:),'rows');

Uages = flip(unique(Ages(Groups==1)),1)

Uages = unique(Ages(Groups==1))
figure;
scatter(Uages,ones(numel(Uages),1),'filled','CData',UcolorsChild)

for u = 1:numel(Uages)
clrs(Ages==Uages(u),:)=repmat(UcolorsChild(u,:),sum(Ages==Uages(u)),1)
end
%}

%%

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