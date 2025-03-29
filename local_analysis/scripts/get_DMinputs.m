function DMinputs = get_DMinputs(SepiUnfilteredMutationVars,CacnesUnfilteredMutationVars,GoodSepiIsolate,GoodCacnesIsolate)
%% Takes a boolean of good samples both species and
% constructs the input for generating a distance matrix
% because it takes too long to run locally

% initialize output
DMinputs = struct;


%% load unfiltered S. epidermidis mutations and resize

load(SepiUnfilteredMutationVars,'SampleNames_unfiltered','maNT_unfiltered','maf_unfiltered','coverage_unfiltered','Quals_unfiltered','indel_counter_unfiltered')
GoodSepiIsolate = ismember(SampleNames_unfiltered,GoodSepiIsolate);

%% DETERMINE MAJOR ALLELES FROM COUNTS
fprintf(1,'Determining major alleles from counts...\n')

maNT = int8(maNT_unfiltered(:,GoodSepiIsolate));
maf = single(maf_unfiltered(:,GoodSepiIsolate));
coverage = coverage_unfiltered(:,GoodSepiIsolate);
Quals = Quals_unfiltered(:,GoodSepiIsolate);
indel_counter=indel_counter_unfiltered(:,:,GoodSepiIsolate);
indel_cov = squeeze(nansum(indel_counter,1));
fraction_indel = indel_cov./(coverage+indel_cov);
Calls = maNT;

% apply first set of filters to Calls_filtered
params.min_qual_for_call = 45;
params.min_maf_for_call = .79;
params.min_cov = 3;
Calls( Quals < params.min_qual_for_call | maf < params.min_maf_for_call | coverage < params.min_cov ) = 0;

diversity = arrayfun(@(x) numel(unique(Calls(x,Calls(x,:)>0))),1:size(Calls,1));
Variable = diversity>1;
CallsVariable = Calls(Variable,:);
QualsVariable = Quals(Variable,:);
OnlyOneMutation = diversity==2;
CallsSingleMutations = Calls(OnlyOneMutation,:);
QualsSingleMutations = Quals(OnlyOneMutation,:);
fNoCallVariable = mean(CallsVariable==0,1);
fNoCallSingleMutation=mean(CallsSingleMutations==0,1);

DMinputs.SepidermidisATCC12228.diversity=diversity;
DMinputs.SepidermidisATCC12228.Variable = Variable;
DMinputs.SepidermidisATCC12228.CallsVariable = CallsVariable;
DMinputs.SepidermidisATCC12228.QualsVariable = QualsVariable;
DMinputs.SepidermidisATCC12228.OnlyOneMutation = OnlyOneMutation;
DMinputs.SepidermidisATCC12228.CallsSingleMutations = CallsSingleMutations;
DMinputs.SepidermidisATCC12228.QualsSingleMutations = QualsSingleMutations;
DMinputs.SepidermidisATCC12228.fNoCallVariable = fNoCallVariable;
DMinputs.SepidermidisATCC12228.fNoCallSingleMutation=fNoCallSingleMutation;
DMinputs.SepidermidisATCC12228.bool=GoodSepiIsolate;
DMinputs.SepidermidisATCC12228.fraction_indel=fraction_indel;

%% As above, load unfiltered C. acnes mutations and resize
load(CacnesUnfilteredMutationVars,'SampleNames_unfiltered','maNT_unfiltered','maf_unfiltered','coverage_unfiltered','Quals_unfiltered','indel_counter_unfiltered')
GoodCacnesIsolate = ismember(SampleNames_unfiltered,GoodCacnesIsolate);

%% DETERMINE MAJOR ALLELES FROM COUNTS
fprintf(1,'Determining major alleles from counts...\n')

maNT = int8(maNT_unfiltered(:,GoodCacnesIsolate));
maf = single(maf_unfiltered(:,GoodCacnesIsolate));
coverage = coverage_unfiltered(:,GoodCacnesIsolate);
Quals = Quals_unfiltered(:,GoodCacnesIsolate);
indel_counter=indel_counter_unfiltered(:,:,GoodCacnesIsolate);
indel_cov = squeeze(nansum(indel_counter,1));
fraction_indel = indel_cov./(coverage+indel_cov);
Calls = maNT;

% apply first set of filters to Calls_filtered
params.min_qual_for_call = 45;
params.min_maf_for_call = .79;
params.min_cov = 3;

Calls( Quals < params.min_qual_for_call | maf < params.min_maf_for_call | coverage < params.min_cov ) = 0;
diversity = arrayfun(@(x) numel(unique(Calls(x,Calls(x,:)>0))),1:size(Calls,1));
Variable = diversity>1;
CallsVariable = Calls(Variable,:);
QualsVariable = Quals(Variable,:);
OnlyOneMutation = diversity==2;
CallsSingleMutations = Calls(OnlyOneMutation,:);
QualsSingleMutations = Quals(OnlyOneMutation,:);
fNoCallVariable = mean(CallsVariable==0,1);
fNoCallSingleMutation=mean(CallsSingleMutations==0,1);
DMinputs.Pacnes_C1.diversity=diversity;
DMinputs.Pacnes_C1.Variable = Variable;
DMinputs.Pacnes_C1.CallsVariable = CallsVariable;
DMinputs.Pacnes_C1.QualsVariable = QualsVariable;
DMinputs.Pacnes_C1.OnlyOneMutation = OnlyOneMutation;
DMinputs.Pacnes_C1.CallsSingleMutations = CallsSingleMutations;
DMinputs.Pacnes_C1.QualsSingleMutations = QualsSingleMutations;
DMinputs.Pacnes_C1.fNoCallVariable = fNoCallVariable;
DMinputs.Pacnes_C1.fNoCallSingleMutation=fNoCallSingleMutation;
DMinputs.Pacnes_C1.bool=GoodSepiIsolate;
DMinputs.Pacnes_C1.fraction_indel=fraction_indel;


end