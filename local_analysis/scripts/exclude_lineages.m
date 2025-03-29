function [Msubject,MsubjectTime,Mtubes]=exclude_lineages(IT,Msubject,MsubjectTime,Mtubes)
%% Removes lineages which incuded in PHLAME classifier, but not clustered because they only contain 2 isolates
% which means you cannot make a phylogeny
% also removes lineages found on researchers

% get number of lineages classified with phlame

CacnesLineageNumbers = Msubject.CacnesLineageNumbers(1,:);
SepiLineageNumbers   = Msubject.SepiLineageNumbers(1,:);

% sget booleans which are true where isolates are found
% this indicuates good clusters because IT only includes good clusters

LineagesToKeepCacnes = arrayfun(@(x) sum(IT.SpeciesString=="cacnes"&IT.ClusterString==x)>0, CacnesLineageNumbers);
LineagesToKeepSepi = arrayfun(@(x) sum(IT.SpeciesString=="sepi"&IT.ClusterString==x)>0, SepiLineageNumbers);

%% Lineages removed because they come from mock communities

SepiLineagesRemove = [6 45];
CacnesLineagesRemove=[43 116];

LineagesToKeepCacnes = LineagesToKeepCacnes&~ismember(CacnesLineageNumbers,CacnesLineagesRemove);
LineagesToKeepSepi = LineagesToKeepSepi&~ismember(SepiLineageNumbers,SepiLineagesRemove);



%% For cacnes first

% lineages for Sepi stay the same because the N=2 clusters were not
% included in S. epi clustering 
% For C. acnes, remove these lineages from abundance info


Msubject=cacnes_lineage_filter(Msubject,LineagesToKeepCacnes);
MsubjectTime=cacnes_lineage_filter(MsubjectTime,LineagesToKeepCacnes);


Mtubes.unfiltered_CacnesLineageAbundance=Mtubes.CacnesLineageAbundance;
Mtubes.unfiltered_CacnesLineageNumbers=Mtubes.CacnesLineageNumbers;
Mtubes.CacnesLineageAbundance=Mtubes.CacnesLineageNumbers(:,LineagesToKeepCacnes);
Mtubes.CacnesLineageNumbers=Mtubes.CacnesLineageNumbers(:,LineagesToKeepCacnes);

%% Then for S. epi
Msubject=sepi_lineage_filter(Msubject,LineagesToKeepSepi);
MsubjectTime=sepi_lineage_filter(MsubjectTime,LineagesToKeepSepi);


Mtubes.unfiltered_SepiLineageAbundance=Mtubes.SepiLineageAbundance;
Mtubes.unfiltered_SepiLineageNumbers=Mtubes.SepiLineageNumbers;
Mtubes.SepiLineageAbundance=Mtubes.SepiLineageNumbers(:,LineagesToKeepSepi);
Mtubes.SepiLineageNumbers=Mtubes.SepiLineageNumbers(:,LineagesToKeepSepi);

%%



function M=cacnes_lineage_filter(M,include)

M.unfiltered_cacnes_clusterlev_concatenated_conservative=M.cacnes_clusterlev_concatenated_conservative;
M.unfiltered_cacnes_clusterlev_concatenated_loose=M.cacnes_clusterlev_concatenated_loose;
M.unfiltered_CacnesLineageNumbers=M.CacnesLineageNumbers;
M.cacnes_clusterlev_concatenated_conservative=M.cacnes_clusterlev_concatenated_conservative(:,include);
M.cacnes_clusterlev_concatenated_loose=M.cacnes_clusterlev_concatenated_loose(:,include);
M.CacnesLineageNumbers=M.CacnesLineageNumbers(:,include);

function M=sepi_lineage_filter(M,include)

M.unfiltered_sepi_clusterlev_concatenated_conservative=M.sepi_clusterlev_concatenated_conservative;
M.unfiltered_sepi_clusterlev_concatenated_loose=M.sepi_clusterlev_concatenated_loose;
M.unfiltered_SepiLineageNumbers=M.SepiLineageNumbers;
M.sepi_clusterlev_concatenated_conservative=M.sepi_clusterlev_concatenated_conservative(:,include);
M.sepi_clusterlev_concatenated_loose=M.sepi_clusterlev_concatenated_loose(:,include);
M.SepiLineageNumbers=M.SepiLineageNumbers(:,include);
