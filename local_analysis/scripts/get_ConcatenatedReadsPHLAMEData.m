function T = get_ConcatenatedReadsPHLAMEData

%% loads in PHLAME data and puts it into matlab structure

sepi_clusterlev_concatenated_conservative =readtable('data/PHLAME_abundances/conservative_Sepi_acera_lineage_level_frequencies_maxpi=0.35_minprob=0.50_hpdfilter=0.1.csv','NumHeaderLines',0);
cacnes_clusterlev_concatenated_conservative =readtable('data/PHLAME_abundances/conservative_Cacnes_acera_lineage_level_frequencies_maxpi=0.35_minprob=0.50_hpdfilter=0.1.csv','NumHeaderLines',0);

sepi_clusterlev_concatenated_loose =readtable('data/PHLAME_abundances/Sepi_acera_lineage_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
cacnes_clusterlev_concatenated_loose =readtable('data/PHLAME_abundances/Cacnes_acera_lineage_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);

sepi_phylolev_concatenated_conservative =readtable('data/PHLAME_abundances/Sepi_acera_phylogroup_level_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
sepi_subphylolev_concatenated_conservative =readtable('data/PHLAME_abundances/Sepi_acera_subphylogroup_level_loose_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
cacnes_phylolev_concatenated_conservative =readtable('data/PHLAME_abundances/Cacnes_acera_phylogroup_level_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);

sepi_phylolev_concatenated_loose =readtable('data/PHLAME_abundances/Sepi_acera_phylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
sepi_subphylolev_concatenated_loose =readtable('data/PHLAME_abundances/Sepi_acera_subphylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
cacnes_phylolev_concatenated_loose =readtable('data/PHLAME_abundances/Cacnes_acera_phylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
% These data were obtained later than the other ones. 
cacnes_subphylolev_concatenated =readtable('data/PHLAME_abundances/Cacnes_acera_subphylogroup_frequencies_2024.csv','NumHeaderLines',0);

%% Parse sample names in PHLAME data to extract TPs
Names = string(sepi_clusterlev_concatenated_conservative.Var1(2:end));
Names=strrep(Names,"ConcatenatedReads_","");
Names=strrep(Names,"TP_","");
Names=strrep(Names,"_all","");

%% Initialize table
T = table;
T.FullNames = Names;
for i = 1:316
    n = strsplit(Names(i),"_");
    if numel(n)==1
        T.SID(i)=n;
    elseif numel(n)==2
        T.SID(i)=n(1);
        T.TP(i)=n(2);
    elseif numel(n)==3
        T.SID(i)=n(1);
        T.TP(i)=n(2);
        T.Plate(i)=n(3);
    end
end

%% Parse clade names

SepiClusterNames = table2array(sepi_clusterlev_concatenated_conservative(1,2:end));
SepiPhyloNames = table2array(sepi_phylolev_concatenated_conservative(1,2:end));
SepiSubhyloNames = string(sepi_subphylolev_concatenated_conservative.Properties.VariableNames(2:end));
CacnesClusterNames = table2array(cacnes_clusterlev_concatenated_conservative(1,2:end));
CacnesPhyloNames = string(cacnes_phylolev_concatenated_conservative.Properties.VariableNames(2:end));
CacnesSubphyloNames = string(cacnes_subphylolev_concatenated.Properties.VariableNames(2:end));

%% Parse abundances

sepi_clusterlev_concatenated_conservative=table2array(sepi_clusterlev_concatenated_conservative(2:end,2:end));
sepi_clusterlev_concatenated_loose=table2array(sepi_clusterlev_concatenated_loose(2:end,2:end));

cacnes_clusterlev_concatenated_conservative=table2array(cacnes_clusterlev_concatenated_conservative(2:end,2:end));
cacnes_clusterlev_concatenated_loose=table2array(cacnes_clusterlev_concatenated_loose(2:end,2:end));

sepi_phylolev_concatenated_conservative=table2array(sepi_phylolev_concatenated_conservative(2:end,2:end));
sepi_subphylolev_concatenated_conservative=table2array(sepi_subphylolev_concatenated_conservative(:,2:end));
cacnes_phylolev_concatenated_conservative=table2array(cacnes_phylolev_concatenated_conservative(:,2:end));

sepi_phylolev_concatenated_loose=table2array(sepi_phylolev_concatenated_loose(2:end,2:end));
sepi_subphylolev_concatenated_loose=table2array(sepi_subphylolev_concatenated_loose(:,2:end));
cacnes_phylolev_concatenated_loose=table2array(cacnes_phylolev_concatenated_loose(:,2:end));
cacnes_subphylolev_concatenated=table2array(cacnes_subphylolev_concatenated(:,2:end));

%% Add abundances to table

T.sepi_clusterlev_concatenated_conservative=sepi_clusterlev_concatenated_conservative;
T.sepi_clusterlev_concatenated_loose=sepi_clusterlev_concatenated_loose;

T.cacnes_clusterlev_concatenated_conservative=cacnes_clusterlev_concatenated_conservative;
T.cacnes_clusterlev_concatenated_loose=cacnes_clusterlev_concatenated_loose;

T.sepi_phylolev_concatenated_conservative=sepi_phylolev_concatenated_conservative;
T.sepi_subphylolev_concatenated_conservative=sepi_subphylolev_concatenated_conservative;
T.cacnes_phylolev_concatenated_conservative=cacnes_phylolev_concatenated_conservative;

T.sepi_phylolev_concatenated_loose=sepi_phylolev_concatenated_loose;
T.sepi_subphylolev_concatenated_loose=sepi_subphylolev_concatenated_loose;
T.cacnes_phylolev_concatenated_loose=cacnes_phylolev_concatenated_loose;
T.cacnes_phylolev_concatenated_loose=cacnes_phylolev_concatenated_loose;
T.cacnes_subphylolev_concatenated=cacnes_subphylolev_concatenated;

%% Add cladenames to table

T.SepiClusterNames=repmat(SepiClusterNames,316,1);
T.SepiPhyloNames=repmat(SepiPhyloNames,316,1);
T.SepiSubhyloNames=repmat(SepiSubhyloNames,316,1);
T.CacnesClusterNames=repmat(CacnesClusterNames,316,1);
T.CacnesPhyloNames=repmat(CacnesPhyloNames,316,1);
T.CacnesSubphyloNames=repmat(CacnesSubphyloNames,316,1);




