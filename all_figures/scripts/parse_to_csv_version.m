function Mnew = parse_to_csv_version(M)

%%

Mnew=table;

Mnew.SID=M.SID;
Mnew.TP=M.TP;
Mnew.location=M.location;
Mnew.Age=M.Age;
Mnew.Sex=M.Sex;

%%

Mnew.CacnesPhylotypeNames=repmat('A B C D E F H K L',size(M,1),1);
Mnew.SepiPhylotypeNames=repmat('1 2 3 4',size(M,1),1);
Mnew.SepiSubPhylotypeNames=repmat('C_1_1 C_1_2_1 C_1_2_2 C_2_1_1 C_2_1_2 C_2_2_1_1 C_2_2_1_2_1_1 C_2_2_1_2_1_2 C_2_2_1_2_1_3 C_2_2_1_2_2_2 C_2_2_1_3_1 C_2_2_1_3_2_1 C_2_2_1_3_2_2_1_1_1 C_2_2_1_3_2_2_2 C_2_2_1_3_2_3 C_2_2_1_3_2_4_1',size(M,1),1);
Mnew.SepiLineageNumbers = repmat(strjoin(string(M.SepiLineageNumbers(1,:)),' '),size(M,1),1);
Mnew.CacnesLineageNumbers = repmat(strjoin(string(M.CacnesLineageNumbers(1,:)),' '),size(M,1),1);
Mnew.SpeciesNamesBracken = repmat(strjoin(string(M.NamesSpecies(1,:)),' '),size(M,1),1);
Mnew.GenusNamesBracken = repmat(strjoin(string(M.NamesGenus(1,:)),' '),size(M,1),1);

%%

Mnew.subject_cutotype=M.subject_cutotype;

%%

R = size(M,1);

Variables = ["BrackenGenus"	"BrackenSpecies"	"CacnesLineageAbundance"	"CacnesPhylotypeAbundance"	"SepiLineageAbundance"	"SepiPhylotypeAbundance"	"subject_cutotype"	"unfiltered_CacnesLineageAbundance"	"unfiltered_CacnesLineageNumbers"	"unfiltered_SepiLineageAbundance"	"unfiltered_SepiLineageNumbers" "sepi_clusterlev_concatenated_conservative"	"sepi_clusterlev_concatenated_loose"	"cacnes_clusterlev_concatenated_conservative"	"cacnes_clusterlev_concatenated_loose"	"sepi_phylolev_concatenated_conservative"	"sepi_phylolev_concatenated_loose"	"cacnes_phylolev_concatenated_loose"	"cacnes_phylolev_concatenated_conservative"	"sepi_subphylolev_concatenated_loose"	"sepi_subphylolev_concatenated_conservative"	"umap_clusters_subjecttime_bracken_species"	"unfiltered_cacnes_clusterlev_concatenated_conservative"	"unfiltered_cacnes_clusterlev_concatenated_loose"	"unfiltered_sepi_clusterlev_concatenated_conservative"	"unfiltered_sepi_clusterlev_concatenated_loose"];

for i = 1:numel(Variables)
    try
    col =M.(Variables(i));
    catch
        continue
    end
    if iscell(col)
        col = vertcat(col{:})
    end
    col_new = arrayfun(@(x) {strjoin(string(col(x,:)),' ')}, 1:R);

    Mnew.(Variables(i))=vertcat(col_new{:});
end


