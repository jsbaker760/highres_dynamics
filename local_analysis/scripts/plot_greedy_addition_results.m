function f=plot_greedy_addition_results(Species)

f=figure;
subplot(2,1,1)
bool = Species.is_anchor==1;
ttl = "ANCHOR";

X=Species.n_clustered_old;
Y=Species.n_clustered_new;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));


scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeColor','black');
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
plot([0 2500],[0 2500],'LineWidth',1,'color','m')

%
subplot(2,1,2)
bool = Species.is_anchor==0;
ttl = "DBSCAN";

X=Species.n_clustered_old;
Y=Species.n_clustered_new;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));


scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeColor','black');
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
xlabel('N clustered beforegreedy addition')
ylabel('N clustered after greedy addition')
title(ttl)

plot([0 2500],[0 2500],'LineWidth',1,'color','m')
pbaspect([1 1 1])