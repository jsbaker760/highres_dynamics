function [anchors, clusters] =  AnchorClust(D,seed,C)
%% ANCHORCLUST
A=(2*C)+1;
D(eye(size(D))==1)=NaN;

anchors = seed;
gaining = true;
while gaining
    idx = ~anchors;
    d2anchors = D(idx,anchors);
    d2closest = min(d2anchors,[],2,'omitnan');
    maxclosest = max(d2closest);

    if maxclosest>=A
        hasmaxclosest = d2closest == maxclosest;
        ia = find(idx);
        ia = ia(hasmaxclosest);
        ia = ia(1);
        anchors(ia)=true;
    else
        gaining = false;
    end
end

clusters = ACradiate_OG(D,anchors, C);

end





