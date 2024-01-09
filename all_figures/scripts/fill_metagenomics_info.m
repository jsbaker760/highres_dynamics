function MetagenomicsSampleTable = fill_metagenomics_info

%% 30 Oct 2023, This version of the script takes the old tube-wise version of metagenomics sample info.'
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


%% split the samplenames for parsing


MetagenomicsSampleNames=arrayfun(@(x) {strsplit(MetagenomicsSampleNames(x),'_')},GoodSamples);
MetagenomicsSampleNames=vertcat(MetagenomicsSampleNames{:});


%% add all information from samplename into table


MetagenomicsSampleTable = table;


MetagenomicsSampleTable.tube =arrayfun(@(x) str2double(x{:}(1:end-1)),cellstr(MetagenomicsSampleNames(:,2)));
MetagenomicsSampleTable.location =arrayfun(@(x) string(x{:}(end)),cellstr(MetagenomicsSampleNames(:,2)));
%
MetagenomicsSampleTable.plate= MetagenomicsSampleNames(:,3);
%
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


%% Biological sex as given in questionairre 


MetagenomicsSampleTable.Sex = arrayfun(@(x) string(unique(SampleInfo.SEX(SampleInfo.SID==x))), MetagenomicsSampleTable.SID);


%% the index corrosponds to the original index (in the list iwth 604 tubes)


MetagenomicsSampleTable.index = GoodSamples;


%% Add  tube data for C. acnes data to table 


% cluster level
clusterlev_tubes = readtable('PHLAME/Cacnes_TUBES_acera_lineage_frequencies_maxpi=0.25_minprob=0.70.csv','NumHeaderLines',0);
clusterlev_tubes = table2array(clusterlev_tubes(2:end,2:end));
clusterlev_tubes=clusterlev_tubes(GoodSamples,:);
Over1 = sum(clusterlev_tubes,2)>1;
Normed = clusterlev_tubes(Over1,:)./sum(clusterlev_tubes(Over1,:),2);
clusterlev_tubes(Over1,:)=Normed;
MetagenomicsSampleTable.CacnesLineageAbundance = clusterlev_tubes;