function DMs=calculate_distance_matrices(DMinputs)
%% calculate pairwise SNP distance matrices for all isolate genomes
% Initialize parallel pools
dbstop if error
distcomp.feature( 'LocalUseMpiexec', false );
restoredefaultpath
availableGPUs = gpuDeviceCount;
pc = parcluster('local');
p=parpool(pc,availableGPUs,'IdleTimeout',Inf);
fprintf('loading DMinputs...\n')

%% hard-coded filters for each of the types of distance matrix to make
% various cutoffs were used, however, only one was used for each species
% the others were used for another project
DMparams=struct;
% turned off
DMparams.fNoCall=[.3 .5 1];
DMparams.FractionN=[.15 .5 1];
ReferenceGenomes = fieldnames(DMinputs);

%% Initialize distance matrix structure
DMs=struct;
fprintf('iterating through reference genomes...\n')
for r = 1:numel(ReferenceGenomes)
    fprintf(['Reference Genome: ' [ReferenceGenomes{r}] '...\n'])

    Inputs = DMinputs.(ReferenceGenomes{r});
    DMs.(ReferenceGenomes{r}) = get_DMs(Inputs,DMparams);
end
end
