function Distances = jaccard(communities)

%% get indices
[NSamples, C]= size(communities);

I =reshape(communities,NSamples,1,C);
J =reshape(communities,1,NSamples,C);

%% Calculate distances
Intersection = sum(I>0&J>0,3);
Union=sum((I>0|J>0),3);

Distances = 1-(Intersection./Union);

%% Remove all-zero comparisons
NullSample = sum(communities,2)==0;
Distances(NullSample|NullSample')=NaN;


%% special case: you are only comparing two things return the scalar value instead of a 2x2
if size(communities,1)==2
    Distances=Distances(2);
end
%%




