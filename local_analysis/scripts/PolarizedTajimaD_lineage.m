function [T] = PolarizedTajimaD_lineage(Mutations,min_columns)

%% first get mean dMRCA

if isempty(Mutations)||size(Mutations,2)<min_columns
    T=NaN;
    N=NaN;
else
    mDMRCA = mean(sum(Mutations,1));
    % then get pairwise distances
    s = size(Mutations,2);
    DM = nan(s);
    for i=1:s
        for j=1:s
            if i~=j
                DM(i,j)=sum(sum([Mutations(:,i) Mutations(:,j)],2)==1);
            end
        end
    end

    % mean pairwise distance
    mpairD=nanmean(DM(:));
    % polarized Tajima's D
    T = mpairD./(2.*mDMRCA);

end
