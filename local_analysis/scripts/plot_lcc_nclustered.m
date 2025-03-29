function f=plot_lcc_nclustered(Species)
%%
f=figure
%%
X=cellfun(@(x) x(1), Species.margins_filter_both);
Y=Species.n_clustered_new;
coefs = cellfun(@(x) x(1), Species.coefs_filter_both);
clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));
%%

subplot(2,1,1)
bool = Species.is_anchor==1;
ttl = "ANCHOR";

% all sub-optimal instances
scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeColor','black');
% if this algorithm has the overall best
if any(bool&Species.instance_used)
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    highest_in_algorithm_used = bool&Species.instance_used;
    hold on
    scatter(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
    scatter(X(highest_in_algorithm_used),Y(highest_in_algorithm_used),'filled','MarkerFaceColor','cyan','MarkerEdgeColor','black');
else
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    hold on
    scatter(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
end
xlabel('cutoff Margin')
ylabel('N clustered')
title(ttl)
pbaspect([1 1 1])
%%
subplot(2,1,2)
bool = Species.is_anchor==0;
ttl = "DBSCAN";

% all sub-optimal instances
scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeColor','black');
% if this algorithm has the overall best
if any(bool&Species.instance_used)
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    highest_in_algorithm_used = bool&Species.instance_used;
    hold on
    scatter(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
    scatter(X(highest_in_algorithm_used),Y(highest_in_algorithm_used),'filled','MarkerFaceColor','cyan','MarkerEdgeColor','black');
else
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    hold on
    scatter(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
end
xlabel('cutoff Margin')
ylabel('N clustered')
title(ttl)
pbaspect([1 1 1])