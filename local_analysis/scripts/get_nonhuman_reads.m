function AllClusteringReads = get_nonhuman_reads(Mtubes,MsubjectTime,RawBrackenDataAllTubes)
%% added to pull out the number of filtered non-human reads used in clustering facial cutotypes for each sample

AllIndividualTubes = string(vertcat(Mtubes.SampleNames{:}));
idxAllIndividualTubes=vertcat(Mtubes.index604{:});

% boolean which is true where the taxon was used in facial cutotype
% clustering and after filtering
taxon_used_in_clustering = ismember(RawBrackenDataAllTubes.SpeciesArray.TaxaNames,MsubjectTime.BrackenNamesClustering(1,:));
M = size(MsubjectTime,1);
AllClusteringReads=zeros(M,1);
for i = 1:M
    index_all_tubes = cell2mat(MsubjectTime.index{i});
    is_sample = ismember(idxAllIndividualTubes,index_all_tubes);
    AllClusteringReads(i) = sum(sum(RawBrackenDataAllTubes.SpeciesArray.new_est_reads(taxon_used_in_clustering,is_sample),2));
end

%%

['mean filtered non-human species-assigned reads/sample=' char(string(round(mean(AllClusteringReads),1,'decimal'))) ', range ' char(string(round(min(AllClusteringReads),1,'decimal'))) '-' char(string(round(max(AllClusteringReads),1,'decimal')))]