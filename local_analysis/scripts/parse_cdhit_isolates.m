function CDHITData=parse_cdhit_isolates(CDHITOutput,AssemblyStats,ncbi_references)
%% SampleNamesCDHIT is the names of all samples which are clustered, isolates as well as reference genomes from NCBI

SampleNamesCDHIT= [string(AssemblyStats.SampleID) ; string(ncbi_references.ncbi_reference_name)];

%% load CDHIT output table

T = readtable(CDHITOutput,'FileType','text','Delimiter',{' ','\t'},'NumHeaderLines',0);

%% Var1 is the number of the homologue cluster (0 indexed)

IsCluster = isnan(T.Var1);

%% Var3 is the samplename which each homologue came from

RowSampleNames       = string(T.Var3);
RowSampleNames = strrep(RowSampleNames,'>','');
RowSampleNames = strrep(RowSampleNames,'...','');

%% scalars used to initialize arrays

C = sum(IsCluster); % number of unique clusters
S = numel(SampleNamesCDHIT); % number of samples
R = numel(IsCluster); % total number of rows in file

%% create I, the number of the homologous cluster (starting at 1, not zeros) for each row in T
% I is 0 where row is not an entry

I = zeros(R,1);
idxC = find(IsCluster);

for r = 1:(numel(idxC)-1)
    I(idxC(r):(idxC(r+1)-1))=r;
end
I(idxC(end):end) = numel(idxC);

I(idxC) = 0;

T.I = I;

%% Create J, the SampleName to which this homologue belongs

% ic gives you the index for each row
[uRowNames,~,ic]=unique(RowSampleNames);

% this gives you the samplename index for row in urows
icdhit_to_sn = arrayfun(@(x) max([0 find(x==SampleNamesCDHIT)]), uRowNames);

% re-indexes such that each entry is index of a samplename in
% SampleNamesCDHIT
J = icdhit_to_sn(ic);


%% create array of size N homologues Clusters x N samples which is true where sample has at least 1 homologue

CDHITClusterSamples = false(C,S);

good_indices = I>0&J>0;
indices_with_homologue= sub2ind([C S],I(good_indices),J(good_indices));
CDHITClusterSamples(indices_with_homologue)=true;

%% add data to structure

CDHITData = struct;
CDHITData.CDHITClusterxSamples = CDHITClusterSamples;
CDHITData.SampleNames = SampleNamesCDHIT;
end

