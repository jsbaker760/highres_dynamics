function clrs = getLinearColors(A,BaseColor,Nbins)

% Diff = max(A)-min(A);
% Range = min(A):(Diff/Nbins):max(A);

Range = 0:(1/Nbins):1;
Range=Range(2:end);
Ucolors = cbrewer2(BaseColor,numel(Range));
% 
% Range(1)=[]
% Ucolors(1,:)=[]

idx = arrayfun(@(x) find(abs(x-Range)==min(abs(x-Range))), A);

clrs = Ucolors(idx,:);
% figure;
% scatter(A,ones(numel(A),1),'filled','CData',clrs)

%
% figure
scatter(1:20,ones(20,1),100,'filled','CData',Ucolors,'Marker','square')
xticks([0 5 10 15 20])
xticklabels(0:.25:1)

end