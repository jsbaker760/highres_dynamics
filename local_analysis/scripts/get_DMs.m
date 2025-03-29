function DMoutputs = get_DMs(DMinputs,DMparams)
%% Calculate distance matrices for both species across a range of parameters
% initialize
dbstop if error
DMoutputs= struct;
[DMoutputs.DMsVariable,DMoutputs.DMsSingle, DMoutputs.excludedSamplesVariable, DMoutputs.excludedSamplesSingle]= deal(cell(numel(DMparams.fNoCall),numel(DMparams.FractionN)));

%% Loop through parameters
for i  = 1:numel(DMparams.fNoCall)
    fprintf(['Excluding bad samples with fNoCall filter ' char(string(i)) ' of ' char(string(numel(DMparams.fNoCall))) '...\n'])

    % find samples with fraction of bad calls
    excludeVariable = DMinputs.fNoCallVariable>DMparams.fNoCall(i);
    excludeSingle = DMinputs.fNoCallSingleMutation>DMparams.fNoCall(i);

    Calls_distVariable = DMinputs.CallsVariable(:,~excludeVariable);
    Calls_distSingle = DMinputs.CallsSingleMutations(:,~excludeSingle);

    FractionNVariable = parallel.pool.Constant(mean(Calls_distVariable==0,2));
    FractionNSingle = parallel.pool.Constant(mean(Calls_distSingle==0,2));

    FractionNcutoffs = parallel.pool.Constant(DMparams.FractionN);

    Calls_distVariable = parallel.pool.Constant(Calls_distVariable);
    Calls_distSingle = parallel.pool.Constant(Calls_distSingle);

    [DMsVariable, DMsSingle]=deal(cell(1,numel(DMparams.FractionN)));
    fprintf(['Running parallel pool ' char(string(i)) ' of ' char(string(numel(DMparams.fNoCall))) '...\n'])
    parfor j = 1:numel(DMparams.FractionN)
        % remove positions based on fraction N
        Calls_distNVariable = Calls_distVariable.Value(FractionNVariable.Value<FractionNcutoffs.Value(j),:);
        Calls_distNSingle = Calls_distSingle.Value(FractionNSingle.Value<FractionNcutoffs.Value(j),:);
        % add to structure
        DMsVariable{1,j}=dm_run(Calls_distNVariable);
        DMsSingle{1,j}=dm_run(Calls_distNSingle);
    end
    DMoutputs.DMsVariable(i,:)=DMsVariable;
    DMoutputs.DMsSingle(i,:)=DMsSingle;
    DMoutputs.excludedSamplesVariable(i,:) = {excludeVariable};
    DMoutputs.excludedSamplesSingle(i,:) = {excludeSingle};
end

end
