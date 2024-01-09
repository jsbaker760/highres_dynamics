function LineageDMRCAsForSupplement = LineageData2supplement(LineageDMRCAsForSupplement)

LineageDMRCAsForSupplement = removevars(LineageDMRCAsForSupplement, ["RTT_clade","RTT_clade_ugenotypes","TreeSNPs","Subjects","UniqueSNPs_genotype"]);
LineageDMRCAsForSupplement = removevars(LineageDMRCAsForSupplement, ["SubjectTreeSNPs_unique_genotypes","SubjectTreeSNPs_isolates"]);

LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates)
LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates=vertcat(LineageDMRCAsForSupplement.UniqueSNPs_subject_isolates{:});
%
LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes)
LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes=vertcat(LineageDMRCAsForSupplement.UniqueSNPs_subject_ugenotypes{:});
%
LineageDMRCAsForSupplement.RTT_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.RTT_subject_isolates)
LineageDMRCAsForSupplement.RTT_subject_isolates=vertcat(LineageDMRCAsForSupplement.RTT_subject_isolates{:});
%
LineageDMRCAsForSupplement.dMRCA_subject_isolates=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.dMRCA_subject_isolates)
LineageDMRCAsForSupplement.dMRCA_subject_isolates=vertcat(LineageDMRCAsForSupplement.dMRCA_subject_isolates{:})
%
LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes=arrayfun(@(x) {strjoin(string(x{:}),' ')}, LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes)
LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes=vertcat(LineageDMRCAsForSupplement.dMRCA_subject_ugenotypes{:});