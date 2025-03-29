function [SampleNames_sepi, SampleNames_cacnes] =  core_genome_filter(CDHITData,AssemblyStats,ncbi_references)
%% get booleans which are true where samplename in CDHIT is a reference genome from NCBI for that species

is_SepidermidisReferences = ismember(CDHITData.SampleNames,ncbi_references.ncbi_reference_name(ncbi_references.name_genus_species=="Staphylococcus epidermidis"));
is_ScapitisReferences = ismember(CDHITData.SampleNames,ncbi_references.ncbi_reference_name(ncbi_references.name_genus_species=="Staphylococcus capitis"));
is_SsaccharolyticusReferences = ismember(CDHITData.SampleNames,ncbi_references.ncbi_reference_name(ncbi_references.name_genus_species=="Staphylococcus saccharolyticus"));
is_CacnesReferences = ismember(CDHITData.SampleNames,ncbi_references.ncbi_reference_name(ncbi_references.name_genus_species=="Cutibacterium acnes"));
is_isolate = ismember(CDHITData.SampleNames,AssemblyStats.SampleID);

%% seperate CDHIT data into data for reference genomes and isolates

CDHIT_Sepidermidis=CDHITData.CDHITClusterxSamples(:,is_SepidermidisReferences);
CDHIT_Scapitis=CDHITData.CDHITClusterxSamples(:,is_ScapitisReferences);
CDHIT_Ssaccharolyticus=CDHITData.CDHITClusterxSamples(:,is_SsaccharolyticusReferences);
CDHIT_Cacnes=CDHITData.CDHITClusterxSamples(:,is_CacnesReferences);

% make sure the CDHIT data for isolates is in the proper index (redundant)
CDHIT_isolates = CDHITData.CDHITClusterxSamples(:,is_isolate);
idx = arrayfun(@(x) find(x==string(AssemblyStats.SampleID)), CDHITData.SampleNames(is_isolate));
CDHIT_isolates=CDHIT_isolates(:,idx);

%% make a figure showing proportion of references for each species which have a given gene 

f=figure;
subplot(2,2,1);
% get f references w gene for S. epidermidis genomes
ReferencePcntSepi = mean(CDHIT_Sepidermidis,2);
histogram(ReferencePcntSepi(ReferencePcntSepi>0),'BinWidth',.02);xlabel('(nonzero) proportion S. epidermidis references with gene');ylabel('Number of genes')
xlim([0 1])

subplot(2,2,2);
% get f references w gene for C. acnes genomes
ReferencePcntCacnes = mean(CDHIT_Cacnes>0,2);
histogram(ReferencePcntCacnes(ReferencePcntCacnes>0),'BinWidth',.02);xlabel('(nonzero) proportion C. acnes references with gene');ylabel('Number of genes')
xlim([0 1])

subplot(2,2,3)
% get f references w gene for S. capitis genomes
ReferencePcntScap = mean(CDHIT_Scapitis>0,2);
histogram(ReferencePcntScap(ReferencePcntScap>0),'BinWidth',.02);xlabel('(nonzero) proportion S. capitis references with gene');ylabel('Number of genes')
xlim([0 1])

subplot(2,2,4)
% get f references w gene for S. saccharolyticus genomes
ReferencePcntSsac = mean(CDHIT_Ssaccharolyticus>0,2);
histogram(ReferencePcntSsac(ReferencePcntSsac>0),'BinWidth',.02);xlabel('(nonzero) proportion S. saccharolyticus references with gene');ylabel('Number of genes')
xlim([0 1])

f.Position=[1009 298 715 680];

%% define for genome for each species 

is_Sepidermidis_core=ReferencePcntSepi==1&ReferencePcntScap==0&ReferencePcntSsac==0;
is_Cacnes_core=ReferencePcntSepi==0&ReferencePcntScap==0&ReferencePcntSsac==0&ReferencePcntCacnes==1;
is_Scapitis_core = ReferencePcntSepi==0&ReferencePcntScap==1&ReferencePcntSsac==0;

%% get proportion of core genome in each isolate 

proportion_Sepidermidis_core_isolates = sum(CDHIT_isolates(is_Sepidermidis_core,:)>0,1)./sum(is_Sepidermidis_core);
proportion_Scapitis_core_isolates = sum(CDHIT_isolates(is_Scapitis_core,:)>0,1)./sum(is_Scapitis_core);
proportion_Cacnes_core_isolates = sum(CDHIT_isolates(is_Cacnes_core,:)>0,1)./sum(is_Cacnes_core);

%% define which isolates are sufficiently good to call S. epidermidis for clustering
 
SepiGoodBool = (proportion_Sepidermidis_core_isolates>.9&proportion_Cacnes_core_isolates<.01&proportion_Scapitis_core_isolates<.01&AssemblyStats.MeanCov'>2);

%% define which isolates are sufficiently good to call C. acnes for clustering 

CacnesGoodBool=proportion_Cacnes_core_isolates>.87&proportion_Sepidermidis_core_isolates<.1;

%% resize booleans to be compatable with SampleNames of size 7443x1 which were used for aligning to reference genomes

SampleNamesAllIsolates=CDHITData.SampleNames(is_isolate);

SampleNames_sepi = SampleNamesAllIsolates(SepiGoodBool);
SampleNames_cacnes = SampleNamesAllIsolates(CacnesGoodBool);



