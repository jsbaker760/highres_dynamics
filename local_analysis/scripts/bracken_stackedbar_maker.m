function f=bracken_stackedbar_maker(Mi,cutoff)

%% plot bars
f=figure;
A = vertcat(Mi.BrackenSpecies{:});
Has = max(A,[],1)>cutoff;
A = A(:,Has);
Names = Mi.NamesSpecies(1,Has);

b = bar(A,'stacked','FaceColor','flat');
cmap = cbrewer2('Spectral',numel(b));
for i = 1:numel(b)
    b(i).CData=cmap(i,:);
end
l = legend;
l.String=Names;
xticks(1:size(Mi))
lbls=arrayfun(@(x) strjoin([Mi.SID(x) string(Mi.TP(x)) Mi.location(x)],''),1:size(Mi,1));
xticklabels(lbls)
