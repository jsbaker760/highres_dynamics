function Msubject=get_species_human_ratio(Msubject)

% load raw bracken data
RawBrackenData=  Brackenfiles2Structure;

% get proportion of reads assigned to human, c acnes, and s. epi
FracHuman=RawBrackenData.SpeciesArray.fraction_total_reads(RawBrackenData.SpeciesArray.TaxaNames=="Homo sapiens",:);
FracSepi=RawBrackenData.SpeciesArray.fraction_total_reads(RawBrackenData.SpeciesArray.TaxaNames=="Staphylococcus epidermidis",:);
FracCacnes=RawBrackenData.SpeciesArray.fraction_total_reads(RawBrackenData.SpeciesArray.TaxaNames=="Cutibacterium acnes",:);

% get ratios (for all tubes)
SepiOverHuman=FracSepi./FracHuman;
CacnesOverHuman=FracCacnes./FracHuman;

S = size(Msubject,1);
[MSepiOverHuman,MCacnesOverHuman]=deal(zeros(S,1));

% iterate through samples and add ratio for each subject

for s=1:S
    index =Msubject.index{s};
    TPs = Msubject.TP{s};
    % average for tubes that have more than one index (Sequenced twice)
    SepiTubes=cellfun(@(x) mean(SepiOverHuman(x)),index);
    CacnesTubes=cellfun(@(x) mean(CacnesOverHuman(x)),index); 
    % average a timepoint
    SepiTPs=arrayfun(@(x) mean(SepiTubes(TPs==x)),unique(TPs));
    CacnesTPs=arrayfun(@(x) mean(CacnesTubes(TPs==x)),unique(TPs));
    % for subject, average across timepoints
    MSepiOverHuman(s)=mean(SepiTPs);
    MCacnesOverHuman(s)=mean(CacnesTPs);

end

% add to subject table
Msubject.SepiOverHuman=MSepiOverHuman;
Msubject.CacnesOverHuman=MCacnesOverHuman;