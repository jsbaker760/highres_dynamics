function [BrackenAbundances,Taxa] = filterBrackenTubes(BrackenDataAll,M,strTaxonLevel)
%% this index is in reference to the original list 604 tubes

GoodSamples = M.index;
N = numel(GoodSamples)
%% load the right raw data

if string(strTaxonLevel)=="Species"
    B = BrackenDataAll.rawdata.SpeciesArray;
elseif string(strTaxonLevel)=="Genus"
    B = BrackenDataAll.rawdata.GenusArray;
end
%% basic filtering samples which are not part of the study

BrackenAbundances = B.fraction_total_reads(:,GoodSamples)';
Taxa = B.TaxaNames;


%% remove empty taxa rows

bad =sum(BrackenAbundances,1)==0;
Taxa(bad)=[];
BrackenAbundances(:,bad)=[];

%% remove control samples

AbundanceSpecies = BrackenAbundances;

SpeciesNames = Taxa;
%%

[~,MostAbundantTaxa_idx] = sort(AbundanceSpecies,2,'descend');
%%

ranked_species = SpeciesNames(MostAbundantTaxa_idx)

%%
AmountAeromonas = startsWith(SpeciesNames,"Aeromonas");
AmountAeromonas=AbundanceSpecies(:,AmountAeromonas)

%% remove human and Aeromonas reads

bad=startsWith(Taxa,"Homo")|startsWith(Taxa,"Aeromonas");
Taxa(bad)=[];
BrackenAbundances(:,bad)=[];

%% renormalzie and remove empty rows again

BrackenAbundances=BrackenAbundances./sum(BrackenAbundances,2);
BrackenAbundances(isnan(BrackenAbundances))=0;
%
%%
Taxa = repmat(Taxa',N,1);