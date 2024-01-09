function Distances = jaccard(communities)



[NSamples, C]= size(communities);

I =reshape(communities,NSamples,1,C);
J =reshape(communities,1,NSamples,C);

%%
Intersection = sum(I>0&J>0,3);
Union=sum((I>0|J>0),3);

Distances = 1-(Intersection./Union);

%%

NullSample = sum(communities,2)==0;
Distances(NullSample|NullSample')=NaN;


% special case: you are only comparing two things, in which case you just
% return the scalar value instead of a 2x2
if size(communities,1)==2
    Distances=Distances(2);
end
%%




