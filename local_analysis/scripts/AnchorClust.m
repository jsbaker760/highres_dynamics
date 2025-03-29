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

clusters = ACradiate(D,anchors, C);

end


function [clusters] = ACradiate(D,anchors, C1)
%% initialize
clusters   = zeros(size(anchors));
Disolates2anchors = (D(~anchors,anchors));

%% number the anchors first

clusters(anchors) = 1:sum(anchors);

%% uniquely closse
Close = Disolates2anchors<=C1;
nClose = sum(Close,2);

points = find(~anchors);

c1 = nClose ==1;
c1p = Close.*((1:size(Close,2)));
clusters(points(c1)) = sum(c1p(c1,:),2);

Close2many = nClose > 1;
clusters(points(Close2many)) = 0;

end


