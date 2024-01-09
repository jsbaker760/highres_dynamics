function T = get_ConcatenatedReadsPHLAMEData



sepi_clusterlev_concatenated_conservative =readtable('PHLAME_filtered/conservative_Sepi_acera_lineage_level_frequencies_maxpi=0.35_minprob=0.50_hpdfilter=0.1.csv','NumHeaderLines',0);
cacnes_clusterlev_concatenated_conservative =readtable('PHLAME_filtered/conservative_Cacnes_acera_lineage_level_frequencies_maxpi=0.35_minprob=0.50_hpdfilter=0.1.csv','NumHeaderLines',0);

sepi_clusterlev_concatenated_loose =readtable('PHLAME_filtered/samples_concatenated/Sepi_acera_lineage_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
cacnes_clusterlev_concatenated_loose =readtable('PHLAME_filtered/samples_concatenated/Cacnes_acera_lineage_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);


sepi_phylolev_concatenated_conservative =readtable('PHLAME_filtered/samples_concatenated/Sepi_acera_phylogroup_level_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
sepi_subphylolev_concatenated_conservative =readtable('PHLAME_filtered/Sepi_acera_subphylogroup_level_loose_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);
cacnes_phylolev_concatenated_conservative =readtable('PHLAME_filtered/samples_concatenated/Cacnes_acera_phylogroup_level_frequencies_maxpi=0.35_minprob=0.50.csv','NumHeaderLines',0);

sepi_phylolev_concatenated_loose =readtable('PHLAME_filtered/samples_concatenated/Sepi_acera_phylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
sepi_subphylolev_concatenated_loose =readtable('PHLAME_filtered/samples_concatenated/Sepi_acera_subphylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
cacnes_phylolev_concatenated_loose =readtable('PHLAME_filtered/samples_concatenated/Cacnes_acera_phylogroup_level/frequencies_maxpi=0.55_minprob=0.50.csv','NumHeaderLines',0);
% 
%%



Names = string(sepi_clusterlev_concatenated_conservative.Var1(2:end));

Names=strrep(Names,"ConcatenatedReads_","")
Names=strrep(Names,"TP_","")
Names=strrep(Names,"_all","");



%%


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

%%

SepiClusterNames = table2array(sepi_clusterlev_concatenated_conservative(1,2:end));
SepiPhyloNames = table2array(sepi_phylolev_concatenated_conservative(1,2:end));
SepiSubhyloNames = string(sepi_subphylolev_concatenated_conservative.Properties.VariableNames(2:end));
CacnesClusterNames = table2array(cacnes_clusterlev_concatenated_conservative(1,2:end));
CacnesPhyloNames = string(cacnes_phylolev_concatenated_conservative.Properties.VariableNames(2:end));
%%
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

%%

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

%%

T.SepiClusterNames=repmat(SepiClusterNames,316,1)
T.SepiPhyloNames=repmat(SepiPhyloNames,316,1)
T.SepiSubhyloNames=repmat(SepiSubhyloNames,316,1)
T.CacnesClusterNames=repmat(CacnesClusterNames,316,1)
T.CacnesPhyloNames=repmat(CacnesPhyloNames,316,1)

%%



