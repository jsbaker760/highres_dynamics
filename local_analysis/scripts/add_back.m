function clusters= add_back(distance_matrix, clusters,maxW)

%% Greedy addition algorithm

clusters(clusters<1)=0;
C = max(clusters);
[closest,odx] = get_closest_o(distance_matrix,clusters);
while numel(odx)>0
    for o=1:numel(odx)
        cidx = odx(o);
        cclust = closest(o);
        cluster_candidate = clusters;
        cluster_candidate(cidx)=cclust;
        row = distance_matrix(cidx,:);
        same=cluster_candidate==cclust;
        clustered=cluster_candidate>0;
        dif = clustered&~same;
        if (all(row(same)<=maxW))
            if max(row(same)) < min(row(dif))
                is_consistent = true;
                d2new = distance_matrix(:,cidx);
                for c = 1:C
                    isc = cluster_candidate==c;
                    if sum(isc)>1 & c~=cclust
                        if any(d2new(isc)<max(distance_matrix(isc',isc),[],2))
                            is_consistent = false;
                            break
                        end
                    end
                end
                if is_consistent
                    clusters = cluster_candidate;
                    [closest,odx] = get_closest_o(distance_matrix,clusters);
                    break
                else
                    clusters(cidx)=-1;
                    odx(o)=[];
                    closest(o)=[];
                    break
                end
            else
                clusters(cidx)=-1;
                odx(o)=[];
                closest(o)=[];
                break
            end
        else
            clusters(cidx)=-1;
            odx(o)=[];
            closest(o)=[];
            break
        end
    end
end
end

function [closest,odx]= get_closest_o(dm,clusters)
clustered = clusters>0;
odx = clusters==0;
clusters = clusters(clustered);
[Dclosest,idxclosest] = min(dm(odx,clustered),[],2);
closest= clusters(idxclosest);
distances_sorted=Dclosest;
odx=find(odx);
if ~issorted(distances_sorted)
    [~,idx]=sort(distances_sorted,'ascend');
    odx = odx(idx);
    closest = closest(idx);
end
end