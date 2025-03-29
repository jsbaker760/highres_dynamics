function f=plot_clusterparams(Species)

%%

f=figure;

%%

subplot(2,2,1)
bool = Species.is_anchor==1;
ttl = "ANCHOR";

X=Species.C1;
Y=Species.n_clustered_old;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));

scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeAlpha',0);
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

xlabel('C1')
ylabel('N clustered')
title(ttl)
pbaspect([1 1 1])


%


subplot(2,2,2)
bool = Species.is_anchor==1;
ttl = "ANCHOR";

X=Species.C1;
Y=Species.n_clustered_new;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));

scatter(X(bool),Y(bool),'filled','CData',clrs(bool,:),'MarkerEdgeAlpha',0);
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

xlabel('C1')
ylabel('N clustered')
title(ttl)
pbaspect([1 1 1])


%


subplot(2,2,3)
bool = Species.is_anchor==0;
ttl = "DBSCAN";

X=Species.epsilon;
Y=Species.minpts;
Z=Species.n_clustered_old;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));


scatter3(X(bool),Y(bool),Z(bool),'filled','CData',clrs(bool,:),'MarkerEdgeAlpha',0);
if any(bool&Species.instance_used)
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    highest_in_algorithm_used = bool&Species.instance_used;
    hold on
    scatter3(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused)-.075,Z(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
    scatter3(X(highest_in_algorithm_used),Y(highest_in_algorithm_used)-.075,Z(highest_in_algorithm_used),'filled','MarkerFaceColor','cyan','MarkerEdgeColor','black');

else
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    hold on
    scatter3(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused)-.075,Z(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
end
% end
xlabel('Epslion')
ylabel('MinPTs')
zlabel('N')
title(ttl)
pbaspect([1 1 1])
%
subplot(2,2,4)
ttl = "DBSCAN";
X=Species.epsilon;
Y=Species.minpts;
Z=Species.n_clustered_new;

coefs = cellfun(@(x) x(1), Species.coefs_filter_both);

clrs = 1-repmat(coefs,1,3);
clrs = clrs./max(clrs(:));


scatter3(X(bool),Y(bool),Z(bool),'filled','CData',clrs(bool,:),'MarkerEdgeAlpha',0);
if any(bool&Species.instance_used)
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    highest_in_algorithm_used = bool&Species.instance_used;
    hold on
    scatter3(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused)-.0,Z(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
    scatter3(X(highest_in_algorithm_used),Y(highest_in_algorithm_used)-.0,Z(highest_in_algorithm_used),'filled','MarkerFaceColor','cyan','MarkerEdgeColor','black');

else
    hold on
    highest_in_algorithm_unused = (bool&Species.instance_best_algo)&~(bool&Species.instance_used);
    hold on
    scatter3(X(highest_in_algorithm_unused),Y(highest_in_algorithm_unused)-.0,Z(highest_in_algorithm_unused),'filled','MarkerFaceColor','magenta','MarkerEdgeColor','black');
end
% end
xlabel('Epslion')
ylabel('MinPTs')
zlabel('N')
title(ttl)
pbaspect([1 1 1])