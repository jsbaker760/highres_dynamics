function f=plot_within_cluster_distance(LD,ManuscriptColors)

%%


L=size(LD,1);
dms=cell(L,1);

for i = 1:L
    ingroup=~LD.outgroup{i};
    S = sum(ingroup);
    muts = LD.calls_analysis{i};
    if isempty(muts)
        dms{i}=zeros(S);
    else
        muts=muts(:,ingroup);
        dm = zeros(S);
        for s = 1:S
            sample_calls = repmat(muts(:,s),1,S);
            dm(s,:)=sum(muts~=sample_calls & muts>0 & sample_calls>0);
        end
        % prevent double-comparison and self-comparison by takng lowe
        % triangle
%         dm(~tril(ones(S),-1)==1)=NaN;
        dms{i}=dm;
    end
end

%%

SepiDistances = dms(LD.SpeciesName=="sepi");
CacnesDistances = dms(LD.SpeciesName=="cacnes");

%% 
SepiDistances=cellfun(@(x) {max(x,[],2,'omitnan')}, SepiDistances);
CacnesDistances=cellfun(@(x) {max(x,[],2,'omitnan')}, CacnesDistances);
%%
SepiDistances=vertcat(SepiDistances{:});
CacnesDistances=vertcat(CacnesDistances{:});
%%
% bw=10;
% f=figure
% histogram(SepiDistances,'BinWidth',bw,'Normalization','probability','FaceColor',ManuscriptColors.SpeciesSepi)
% hold on
% histogram(CacnesDistances,'BinWidth',bw,'Normalization','probability','FaceColor',ManuscriptColors.SpeciesCacnes)
% xlim([0 100])

%

[sepibins,edges]= histcounts(SepiDistances,0:bw:100);
[cacnesbins,edges]= histcounts(CacnesDistances,0:bw:100);

f=figure;
b=bar(([sepibins ; cacnesbins])','stacked','FaceColor','flat')
b(1).CData = [0 0 0];
b(2).CData = [1 1 1];

xticks(1:20);
xticklabels(edges(2:end))

xlabel('Maximum within-cluster distance(lineage co-assemblies)')
ylabel('N isolates')
%

text(0,1400,['maxCacnes=' char(string(max(CacnesDistances)))])
text(0,1200,['maxSepi=' char(string(max(SepiDistances)))])

%%




