function [N,W,C] = getNWC(DM,clusters)
%%
[R,S] = size(clusters);
[W,N,C] = deal(zeros(R,1));

for r = 1:R
    row = clusters(r,:);
    W(r) = max(arrayfun(@(x) max(max(DM(row==x,row==x))),1:max(row)));
    T = tabulate(row(row>0));
    T = T(:,2);
    N(r) = sum(T(T>3));
end
%%
end