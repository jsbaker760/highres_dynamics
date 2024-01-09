function [closest,odx]= get_closest_o(dm,clusters)

%%
clustered = clusters>0;
odx = clusters==0;
clusters = clusters(clustered);
%%
[Dclosest,idxclosest] = min(dm(odx,clustered),[],2);
%%
closest= clusters(idxclosest);
%%
distances_sorted=Dclosest;
odx=find(odx);
%%
if ~issorted(distances_sorted)
    [~,idx]=sort(distances_sorted,'ascend');
    odx = odx(idx);
    closest = closest(idx);
end
