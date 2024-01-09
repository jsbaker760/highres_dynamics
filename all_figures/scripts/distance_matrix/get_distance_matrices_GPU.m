%% pxNsample for each of two reference genomes for sampels considered good
%%
dbstop if error
%%
distcomp.feature( 'LocalUseMpiexec', false );
restoredefaultpath
availableGPUs = gpuDeviceCount;
pc = parcluster('local');
p=parpool(pc,availableGPUs,'IdleTimeout',Inf);

%%
fprintf('loading DMinputs...\n')

%%
 load DMinputs.mat
%% hard-coded filters for each of the types of distance matrix to make;
DMparams=struct;
% turned off
DMparams.fNoCall=[.3 .5 1];
DMparams.FractionN=[.15 .5 1];
%%
ReferenceGenomes = fieldnames(DMinputs);
%%
DMs=struct;
%%
fprintf('iterating through reference genomes...\n')
for r = 1:numel(ReferenceGenomes)
fprintf(['Reference Genome: ' [ReferenceGenomes{r}] '...\n'])

    Inputs = DMinputs.(ReferenceGenomes{r});
    DMs.(ReferenceGenomes{r}) = get_DMs(Inputs,DMparams);
end
%%
save('DMs.mat','DMs','DMparams','-v7.3')