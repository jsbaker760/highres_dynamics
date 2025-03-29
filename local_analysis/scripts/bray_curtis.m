function Distances = bray_curtis(communities)

if numel(size(communities))>2
  error('input is supposed to be 2-dimentional (Nsamples x Number of Taxa)')
end


[NSamples, C]= size(communities);


%%

I =repmat(reshape(communities,NSamples,1,C),1,NSamples,1);
J =repmat(reshape(communities,1,NSamples,C),NSamples,1,1);

%%
Both = I>0&J>0;
Lesser = cat(4,I,J);
Lesser = min(Lesser,[],4);
Lesser(~Both)=0;
Cij =sum(Lesser,3);

%%

Si = sum(I ,3);
Sj =sum(J,3);


Distances = 1 - ((2*Cij)./(Si+Sj));

%%
NullSample = sum(communities,2)==0;
Distances(NullSample|NullSample')=NaN;
%%
% special case: you are only comparing two things, in which case you just
% return the scalar value instead of a 2x2
% if all(size(Distances)==[2 2])
%     Distances=Distances(2);
% end

end