function MetagenomicsSampleTable = fill_metagenomics_info_filtered
%% load C. acnes abundance information and add to the table



MetagenomicsSampleNames = string(table2array(readtable('metagenomics_samplenames_tubes.csv','Delimiter',',')));
SampleInfo=readtable('subject_sampling_metadata.csv','Delimiter',',','VariableNamingRule','preserve');




%% turn dates into timepoints


dates=datetime(SampleInfo.SamplingDate);
TP = dates_to_timepoints(dates);

SampleInfo.TP=TP;
% dates=string(SampleInfo.SamplingDate);
% SampleInfo.TP = zeros(K,1);
% SampleInfo.TP(dates=="2018-05-29"|dates=="2018-05-30")=1;
% SampleInfo.TP(dates=="2018-10-25")=2;
% SampleInfo.TP(dates=="2018-12-12"|dates=="2018-12-13")=3;
% SampleInfo.TP(dates=="2019-06-04")=4;
% SampleInfo.TP(dates=="2019-10-24")=5;
% SampleInfo.TP(dates=="2022-12-02")=6;





%% remove samples which are irrelevant 




Irrelevent = contains(MetagenomicsSampleNames,"MOCK")   | ...
             contains(MetagenomicsSampleNames,"NEG")    | ...
             contains(MetagenomicsSampleNames,"BLANK")  | ...
             contains(MetagenomicsSampleNames,"C2C2")   | ...
             contains(MetagenomicsSampleNames,"ISO") ;


GoodSamples = find(~Irrelevent);


G=numel(GoodSamples);


%% split the samplenames for parsing



MetagenomicsSampleNamesFull = MetagenomicsSampleNames;
MetagenomicsSampleNames=arrayfun(@(x) {strsplit(MetagenomicsSampleNames(x),'_')},GoodSamples);
MetagenomicsSampleNames=vertcat(MetagenomicsSampleNames{:});



%% add all information from samplename into table



MetagenomicsSampleTable = table;
% tube is the number actually written on the tube
MetagenomicsSampleTable.tube =arrayfun(@(x) str2double(x{:}(1:end-1)),cellstr(MetagenomicsSampleNames(:,2)));
% Forehead, Nose, cheeK, Chin
MetagenomicsSampleTable.location =arrayfun(@(x) string(x{:}(end)),cellstr(MetagenomicsSampleNames(:,2)));
% Name of 96 well plate used in sequencing (MG1 through MG7)
MetagenomicsSampleTable.plate= MetagenomicsSampleNames(:,3);
% Well number, A1=1, B1=13, H12=96
MetagenomicsSampleTable.well = MetagenomicsSampleNames(:,4);



%% add SID and timepoint to rows



SampleIDused = SampleInfo.SampleNameCode;
MetagenomicsSampleTable.SID = strings(size(MetagenomicsSampleTable.tube));
MetagenomicsSampleTable.TP = zeros(size(MetagenomicsSampleTable.tube));

for m =1:numel(MetagenomicsSampleTable.SID)
    tube = MetagenomicsSampleTable.tube(m);
    idx = tube==MetagenomicsSampleTable.tube;
    MetagenomicsSampleTable.SID(idx) = SampleInfo.SID(SampleIDused==tube);
    MetagenomicsSampleTable.TP(idx) = SampleInfo.TP(SampleIDused==tube);
    MetagenomicsSampleTable.Age(idx) = years(SampleInfo.SamplingDate(SampleIDused==tube)-SampleInfo.DOB(SampleIDused==tube));
end



%% add sex



MetagenomicsSampleTable.Sex = arrayfun(@(x) string(unique(SampleInfo.SEX(SampleInfo.SID==x))), MetagenomicsSampleTable.SID);



%% the index corrosponds to the original index (in the list iwth 604 tubes)



MetagenomicsSampleTable.index = GoodSamples;



%% Add  tube data for C. acnes data to table 


% Lineage-level abundances
clusterlev_tubes = readtable('PHLAME_filtered/samples_tubewise/Cacnes_acera_lineage_level_tubewise_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
clusterlev_tubes = table2array(clusterlev_tubes(2:end,2:end));
clusterlev_tubes=clusterlev_tubes(GoodSamples,:);
% if abundance is more than one, normalize it
Over1 = sum(clusterlev_tubes,2)>1;
Normed = clusterlev_tubes(Over1,:)./sum(clusterlev_tubes(Over1,:),2);
clusterlev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.CacnesLineageAbundance = clusterlev_tubes;

% phylogroup level
phylolev_tubes = readtable('PHLAME_filtered/samples_tubewise/Cacnes_acera_phylogroup_level_tubewise_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
PhylotypeNames = string(phylolev_tubes.Properties.VariableNames(2:end));
[~,idx]=sort(PhylotypeNames);
phylolev_tubes = table2array(phylolev_tubes(:,2:end));
phylolev_tubes = phylolev_tubes(:,idx);
phylolev_tubes=phylolev_tubes(GoodSamples,:);
Over1 = sum(phylolev_tubes,2)>1;
Normed = phylolev_tubes(Over1,:)./sum(phylolev_tubes(Over1,:),2);
phylolev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.CacnesPhylotypeAbundance = phylolev_tubes;
MetagenomicsSampleTable.CacnesPhylotypeNames = repmat(PhylotypeNames,G,1);

% % add coverage
tubecov = readtable('Cacnes_coverage.txt');
MetagenomicsSampleTable.CacnesMGCoverage = tubecov.Var1(GoodSamples);




%% Add S. epi data to table



 % lineage level
clusterlev_tubes = readtable('PHLAME_filtered/samples_tubewise/Sepi_acera_lineage_level_tubewise_frequencies_maxpi=0.35_minprob=0.50.csv');
clusterlev_tubes = table2array(clusterlev_tubes(2:end,2:end));
clusterlev_tubes=clusterlev_tubes(GoodSamples,:);
Over1 = sum(clusterlev_tubes,2)>1;
Normed = clusterlev_tubes(Over1,:)./sum(clusterlev_tubes(Over1,:),2);
clusterlev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.SepiLineageAbundance = clusterlev_tubes;

%. phylogroup level
phylolev_tubes = readtable('PHLAME_filtered/samples_tubewise/Sepi_acera_phylogroup_level_tubewise_frequencies_maxpi=0.35_minprob=0.50.csv');
phylolev_tubes = table2array(phylolev_tubes(2:end,2:end));
phylolev_tubes=phylolev_tubes(GoodSamples,:);
Over1 = sum(phylolev_tubes,2)>1;
Normed = phylolev_tubes(Over1,:)./sum(phylolev_tubes(Over1,:),2);
phylolev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.SepiPhylotypeAbundance = phylolev_tubes;

%. subsubphylogroup level
subphylolev_tubes = readtable('PHLAME_filtered/samples_tubewise/Sepi_acera_subphylogroup_level_tubewise_frequencies_maxpi=0.35_minprob=0.50.csv');
subphylo_names=string(subphylolev_tubes.Properties.VariableNames(2:end));
subphylolev_tubes = table2array(subphylolev_tubes(:,2:end));
subphylolev_tubes=subphylolev_tubes(GoodSamples,:);
Over1 = sum(subphylolev_tubes,2)>1;
Normed = subphylolev_tubes(Over1,:)./sum(subphylolev_tubes(Over1,:),2);
subphylolev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.SepiSubphylotypeAbundance = subphylolev_tubes;

% genomic coverage
tubecov = readtable('Sepi_coverage.txt');
MetagenomicsSampleTable.SepiMGCoverage = tubecov.Var1(GoodSamples);
% info about which lineages belong to which phylotypes



%% add lineage numbers and phylotypes of each lineage
% add phylotypes of lineages



lineages_to_phylotypes = readtable('PHLAME_filtered/clade_IDs/Cacnes_acera_lineage2phylogroup.txt');
MetagenomicsSampleTable.CacnesLineageNumbers = repmat(lineages_to_phylotypes.Var2',G,1);
MetagenomicsSampleTable.CacnesLineagePhylotypes = repmat(string(lineages_to_phylotypes.Var1)',G,1);

lineages_to_phylotypes = readtable('PHLAME_filtered/clade_IDs/Sepi_acera_lineage2phylogroup.txt');
MetagenomicsSampleTable.SepiLineageNumbers = repmat(lineages_to_phylotypes.Var1',G,1);
MetagenomicsSampleTable.SepiLineagePhylotypes = repmat(string(lineages_to_phylotypes.Var2)',G,1);

MetagenomicsSampleTable.SepiSubphylotypeNames = repmat(subphylo_names,G,1);

%% add original (full) samplenames to each row for posterity



MetagenomicsSampleTable.SampleNames = MetagenomicsSampleNamesFull(GoodSamples);



%% load in raw bracken data from metagenomics samples



BrackenDataAll.rawdata =  Brackenfiles2Structure;



%% add the reads assigned at the genus and species level 



GenusReads = (BrackenDataAll.rawdata.GenusArray.new_est_reads)';
SpeciesReads = (BrackenDataAll.rawdata.SpeciesArray.new_est_reads)';

MetagenomicsSampleTable.BrackenReadsAssignedGenus = GenusReads(GoodSamples,:);
MetagenomicsSampleTable.BrackenReadsAssignedSpecies = SpeciesReads(GoodSamples,:);



%% add species and genus-level abundances from bracken 



[MetagenomicsSampleTable.BrackenSpecies, MetagenomicsSampleTable.NamesSpecies]=filterBrackenTubes(BrackenDataAll,MetagenomicsSampleTable,'Species');
[MetagenomicsSampleTable.BrackenGenus, MetagenomicsSampleTable.NamesGenus]=filterBrackenTubes(BrackenDataAll,MetagenomicsSampleTable,'Genus');


