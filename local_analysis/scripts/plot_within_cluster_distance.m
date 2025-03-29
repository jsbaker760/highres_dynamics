function f=plot_within_cluster_distance(LD)

%% iterate through lineages and get inter-genome intra-lineage distances

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
        dms{i}=dm;
    end
end
% seperate by species
SepiDistances = dms(LD.SpeciesName=="sepi");
CacnesDistances = dms(LD.SpeciesName=="cacnes");
% get longest distance
SepiDistances=cellfun(@(x) {max(x,[],2,'omitnan')}, SepiDistances);
CacnesDistances=cellfun(@(x) {max(x,[],2,'omitnan')}, CacnesDistances);
% put them together for histogram
SepiDistances=vertcat(SepiDistances{:});
CacnesDistances=vertcat(CacnesDistances{:});

%% make the plot
bw=10;

f=figure;

subplot(2,1,1)
histogram(CacnesDistances,'BinWidth',bw,'FaceColor',[1 1 1])
xlim([0 100])
ylim([0 800])
xticks([0:20:100])
yticks([0:200:800])
xlabel('Maximum within-cluster distance(lineage co-assemblies)')
ylabel('N isolates')
CacnesMax=char(string(max(CacnesDistances)))
text(80,600,['Max=' CacnesMax])

subplot(2,1,2)
histogram(SepiDistances,'BinWidth',bw,'FaceColor',[0 0 0])
xlim([0 100])
ylim([0 800])
xticks([0:20:100])
yticks([0:200:800])
xlabel('Maximum within-cluster distance(lineage co-assemblies)')
ylabel('N isolates')
SepiMax=char(string(max(SepiDistances)));
text(80,600,['Max= SepiMax' ])
%%

['Isolates from the same lineage are separated by ≤ ' CacnesMax ' and ≤ ' SepiMax ' SNPs across the whole genome for C. acnes and S. epidermidis, respectively']