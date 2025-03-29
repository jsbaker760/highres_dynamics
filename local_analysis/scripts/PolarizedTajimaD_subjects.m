function [T,N] = PolarizedTajimaD_subjects(Mutations,min_columns)

%% first get mean dMRCA

N=arrayfun(@(x) size(x{:},2), Mutations);
mDMRCA = arrayfun(@(x) mean(sum(x{:},1)), Mutations);

% then get pairwise distances
DMs = cell(size(Mutations));
for m = 1:numel(Mutations)
    clademuts = Mutations{m};
    s = size(clademuts,2);
    dm = nan(s);
    for i=1:s
        for j=1:s
            if i~=j
                dm(i,j)=sum(sum([clademuts(:,i) clademuts(:,j)],2)==1);
            end
        end
    end
    DMs{m}=dm;
end


% mean pairwise distance
mpairD=arrayfun(@(x) mean(x{:}(:),'omitnan'), DMs);

% Polarized Tajima's D

T = mpairD./(2.*mDMRCA);
T(N<min_columns)=NaN;

histogram(T,'BinWidth',.1,'Normalization','probability');
xlim([0 1])
ylim([0 .2])
xlabel('Polarized Tajima D')
ylabel('Proportion subject lineages')
