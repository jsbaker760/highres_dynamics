function UnfilteredMutationVars = resize_candidate_mutation_table(candidate_mutation_table)
% initialize
UnfilteredMutationVars=struct;

% From case step:
fprintf(1,'Loading data...\n')

% get sample names
load(candidate_mutation_table,'SampleNames')
UnfilteredMutationVars.SampleNames_unfiltered = SampleNames';

% counts matrix
load(candidate_mutation_table,'counts')
counts = uint16(counts);
UnfilteredMutationVars.coverage_forward_strand_unfiltered=double(squeeze(sum(counts(1:4,:,:))));
UnfilteredMutationVars.coverage_reverse_strand_unfiltered=double(squeeze(sum(counts(5:8,:,:))));
UnfilteredMutationVars.coverage_unfiltered=coverage_forward_strand_unfiltered+coverage_reverse_strand_unfiltered;
UnfilteredMutationVars.counts_unfiltered = counts;
clear counts

% Quality scores
load(candidate_mutation_table,'Quals')
UnfilteredMutationVars.Quals_unfiltered = int16(-1*Quals);
clear Quals

% positions and other smaller variables
load(candidate_mutation_table,'p','indel_counter')
UnfilteredMutationVars.p_all = p;
UnfilteredMutationVars.LongNames_unfiltered=SampleNames; % original sample names
UnfilteredMutationVars.indel_counter_unfiltered=uint16(indel_counter);

clear p
clear SampleNames
clear indel_counter


% Calculate coverage, call major alleles, filter calls...
fprintf(1,'Calculating coverage from counts...\n')
% DETERMINE MAJOR ALLELES FROM COUNTS
fprintf(1,'Determining major alleles from counts...\n')
[UnfilteredMutationVars.maf_unfiltered, UnfilteredMutationVars.maNT_unfiltered, ~, ~] = div_major_allele_freq(counts_unfiltered);
