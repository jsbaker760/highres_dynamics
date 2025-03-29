function AllIntervalData=get_intervals_oh_et_al(Metadata,LocationsToPlot,CacnesAbundances,SepiAbundances,CacnesSampleNames,SepiSampleNames)

%% initialize indices
UniqueSubjects = unique(Metadata.SubjectNumber); % subject names are parsed into numbers 
L = numel(LocationsToPlot);% number of locations to plot data and therefore subpanels
S = numel(UniqueSubjects);%

%% get temporal data
AllIntervalData = cell(L*S,1);

r = 0;
for l=1:L
    for s = 1:S
        T = table;
        % pull out metadata for all samples from this subject and location
        SubjectLocationMetadata = Metadata(Metadata.SubjectNumber==s&startsWith(Metadata.Site_Symmetry,LocationsToPlot(l)),:);
        if size(SubjectLocationMetadata,1)<2
            continue
        end
        % sort it by timepoint (initial sampling first)
        [~,sortidx]=sort(SubjectLocationMetadata.DateCollected,'ascend');
        SubjectLocationMetadata=SubjectLocationMetadata(sortidx,:);
        idxcacnes = arrayfun(@(x) {find(x==CacnesSampleNames)}, SubjectLocationMetadata.SAMPLE);
        idxsepi = arrayfun(@(x) {find(x==SepiSampleNames)}, SubjectLocationMetadata.SAMPLE);
        missing_cacnes = cellfun(@isempty,idxcacnes);
        missing_sepi = cellfun(@isempty,idxcacnes);
        if any(missing_cacnes|missing_sepi)
            has_both =~missing_sepi&~missing_cacnes;
            if sum(has_both)<2
                continue
            else
                SubjectLocationMetadata = SubjectLocationMetadata(has_both,:); % when one is missing, the other is also always missing
                idxcacnes=idxcacnes(has_both);
                idxsepi=idxsepi(has_both);
            end
        end
        r=r+1;
        idxcacnes=vertcat(idxcacnes{:});
        idxsepi=vertcat(idxsepi{:});
        R = size(SubjectLocationMetadata,1);
        % error out if there are multiple timepoint samples
        if numel(unique(SubjectLocationMetadata.DateCollected))~=R
            error('should only have one entry per subject per location per timepoint')
        end
        % index
        % amount lost
        T.AbundanceLostCacnes = {get_abundance_lost(CacnesAbundances(:,idxcacnes))};
        T.AbundanceLostSepi = {get_abundance_lost(SepiAbundances(:,idxsepi))};
        % bray curtis dissimilarity
        BrayCurtisCacnes = bray_curtis(CacnesAbundances(:,idxcacnes)');
        BrayCurtisSepi = bray_curtis(SepiAbundances(:,idxsepi)');
        T.BrayCurtisCacnes={BrayCurtisCacnes(1,:)};
        T.BrayCurtisSepi={BrayCurtisSepi(1,:)};
        T.IntervalTimes = {years(diff(SubjectLocationMetadata.DateCollected)')};
        T.CumulativeTime={years(cumsum([0 diff(SubjectLocationMetadata.DateCollected)']))};
        T.SubjectNumber = s;
        T.Location = LocationsToPlot(l);
        AllIntervalData{r}=T;

    end
end
AllIntervalData = vertcat(AllIntervalData{1:r});


function Alost= get_abundance_lost(abundance)
    A1 = abundance(:,1:R-1);% starting abundance
    A2 = abundance(:,2:R); % ending abundance
    is_lost = A1>0&A2==0; % true where lost
    A1(~is_lost)=0; % remove non-lost abundances 
    Alost = sum(A1);% add lost abundanced together

    end

end


