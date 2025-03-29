function [clades] = identify_mutations_lineage_assemblies(clustersdir,cladenames,params)
%% INPUTS
% clustersdir is a directory with data for each clade in the format
% "clade_" followed by the number.
% each dir must have: candidate_mutation_table.mat and genome.fasta
% params: 11 different parameters to call mutations
% C number(s) of clusters to look at

C = numel(cladenames);

workingdir=char(pwd);
%% add DNApars to folder where phylogenies are going to be added
if ~exist(clustersdir,'dir')
    mkdir(clustersdir)
    copyfile('dnapars', [clustersdir '/dnapars'])
    copyfile('dnapars.app/', [clustersdir '/dnapars.app/'])
end

treedir = ([workingdir '/' clustersdir]);


%% initialize cell arrays for DMRCA analysis
fprintf('intializing clade structure\n')
clades = struct;
[clades.anc_nti,           ...
    clades.goodpos,            ...
    clades.outgroup,           ...
    clades.calls_analysis,     ...
    clades.samplenames,        ...
    clades.p,                  ...
    clades.fraction_ambiguous, ...
    clades.has_mutation_ingroup_goodpos_no_reco, ...
    clades.annotations,        ...
    clades.tree_snps,...
    clades.outgroup_calls,     ...
    clades.included_after_strict_filtering, ...
    clades.maf,...
    clades.outgroup_cov]       = deal(cell(max(C),1));

clades.cladenumber = zeros(C,1);

%% iterate through clusters
fprintf('starting for-loop\n')

for k = 1:C
    %% get clade name and parse it to get the number


    CladeName = cladenames(k);
    CladeNumer = strsplit(CladeName,"_");
    CladeNumer=str2double(CladeNumer(3));

    %% get the reference genome directory for this clade


    cladedir = strjoin([treedir CladeName],'/');
    REFGENOMEFOLDER = char(strjoin([params.REFGENOMEDIRassemblies CladeName],'/'));


    %% load in reference genome


    [ChrStarts, ~, ~, ~] = genomestats(REFGENOMEFOLDER);


    %% load in candidate_mutation_table.mat


    load(strjoin([REFGENOMEFOLDER "candidate_mutation_table.mat"],'/'),'Quals','SampleNames','counts','in_outgroup','p');


    %% skip this clade if there's no mutations, but still record samplenames and clade number


    in_outgroup = in_outgroup==1;
    SampleNames_ingroup = SampleNames(~in_outgroup);
    clades.cladenumber(k) = CladeNumer;
    clades.samplenames{k} = SampleNames;
    clades.outgroup{k}=in_outgroup;
    if isempty(p)
        continue
    end


    %% inverse quals


    Quals = -1.*Quals;


    %% read in reference genome


    chrpos = p2chrpos(p,ChrStarts);
    refnt_all = extract_outgroup_mutation_positions(REFGENOMEFOLDER, chrpos);
    [~,refnti]=ismember(refnt_all,params.NTs);


    %% make a directory to put this cluster's phylogeny


    % Make a directory for this cluster

    mkdir(cladedir)
    copyfile('dnapars', strjoin([cladedir "dnapars"],'/'))
    copyfile('dnapars.app', strjoin([cladedir "dnapars.app"],'/'))


    %% Nsamples is both ingroup and outgroup


    Nsample=numel(SampleNames);


    %% tag outgroup samplenames for phylogeny


    % Add outgroup tag to a new version of SampleNames
    SampleNames_tagoutgroup = string(SampleNames);
    SampleNames_tagoutgroup(in_outgroup)=arrayfun(@(x) strjoin([ "OUTGROUP" x],'_'),SampleNames(in_outgroup));


    %% get name for phylogeny based on subjects whose isolates are in clade


    subjects_cluster_ingroup = char(SampleNames(~in_outgroup));
    subjects_cluster_ingroup = unique(string(subjects_cluster_ingroup(:,1:3)))';
    subjects_cluster_ingroup = join(subjects_cluster_ingroup,'_');
    subjects_cluster_ingroup_list=subjects_cluster_ingroup;
    subjects_cluster_ingroup = strjoin([char(CladeName) '_subjects_' subjects_cluster_ingroup],'');
    fprintf(1,[char(CladeName) '-subjects_~' subjects_cluster_ingroup{1} '\n'])


    %% create maf and maNT from counts


    [maf, maNT, ~, ~] = div_major_allele_freq(counts);

    clades.maf{k}=maf;

    %% initialize Calls and coverage


    Calls = maNT;
    coverage = squeeze(sum(counts,1));


    %% Filter Calls for SNPs, setting low quality calls to zero


    Calls( Quals     <          params.min_qual_for_call     |...
        maf       <          params.min_maf_for_call      |...
        coverage  <          params.min_cov_for_call      ) = 0;


    %% get N ingroup samples and ingroup calls to remove non-diverse positions


    Nsamples_ingroup = sum(~in_outgroup);
    Calls_ingroup = Calls(:,~in_outgroup);


    %% Find nonvariable positions


    % true where each ingroup sample has the same nucleotide
    nonvariablep=(sum(Calls_ingroup==1,2)==Nsamples_ingroup | sum(Calls_ingroup==2,2)==Nsamples_ingroup | sum(Calls_ingroup==3,2)==Nsamples_ingroup | sum(Calls_ingroup==4,2)==Nsamples_ingroup | sum(Calls_ingroup==0,2)==Nsamples_ingroup);


    %% remove nonvariable positions from data
    fprintf(1,['removing ' num2str(sum(nonvariablep)) 'non-variable positions \n']);


    Calls_ingroup = Calls_ingroup(~nonvariablep,:);
    p=p(~nonvariablep);
    Calls=Calls(~nonvariablep,:);
    counts=counts(:,~nonvariablep,:);
    Quals=Quals(~nonvariablep,:);
    maf=maf(~nonvariablep,:);
    coverage=coverage(~nonvariablep,:);
    maNT=maNT(~nonvariablep,:);
    refnti = refnti(~nonvariablep);
    chrpos = chrpos(~nonvariablep,:);


    %% Define ancestor
    fprintf(1,'Determining outgroup nucleotides...\n');


    anc_nti = nan(size(p));
    % FIRST Define ancestor based on outgroup
    outgroup_calls = Calls(:,in_outgroup);
    outgroup_calls( outgroup_calls==0 ) = nan;  %mode takes the lowest value when there is a tie, setting to nan avoids this
    unambig_outgroup = arrayfun(@(x) numel(unique(outgroup_calls(x,~isnan(outgroup_calls(x,:))))) , 1:numel(p))'==1;
    n_outgroups_with_call= sum(outgroup_calls>0,2);
    use_outgroup_positions=unambig_outgroup&n_outgroups_with_call>1;
    anc_nti(use_outgroup_positions) = mode(outgroup_calls(use_outgroup_positions,:),2);
    fprintf(1,[num2str(sum(isnan(anc_nti))) ' positions without outgroup call\n']);

    % SECOND Replace with reference if no value exists from outgroup
    anc_nti(isnan(anc_nti)) = refnti(isnan(anc_nti));
    fprintf(1,[num2str(sum(isnan(anc_nti))) ' positions without reference genome call\n']);

    % THIRD Use mode of ingroup calls for locations without outgroup or reference
    locs_without_ancestor = sum(Calls(:,~in_outgroup)==repmat(anc_nti,1,Nsamples_ingroup),2)==0;
    ingroup_anc = Calls(locs_without_ancestor,~in_outgroup);
    ingroup_anc(ingroup_anc==0) = nan;  %mode takes the lowest value when there is a tie, setting to nan avoids this
    anc_nti(locs_without_ancestor) = mode(ingroup_anc,2);


    %% remove positions from calls where too ambiguous or not enough coverage
    fprintf(1,'Filtering positions and identifying fixed mutations...\n');


    fraction_ambiguous_ingroup=mean(Calls_ingroup==0,2);
    clades.fraction_ambiguous{k} = fraction_ambiguous_ingroup;
    Calls( fraction_ambiguous_ingroup > params.max_fraction_ambigious_samples,:) = 0;


    %% Define MutQual and fixedmutation with outgroup isolates for tree analysis


    [MutQual_all, ~] = ana_mutation_quality( Calls,Quals );
    fixedmutation_tree = ( (Calls~=repmat(anc_nti,1,Nsample)) & Calls~=0 ); % fixed mutations where call is different from the ancestor and where the call is not 0
    fixedmutation_tree( MutQual_all < 1,: ) = 0; % ...and not where there is only one type of non-N call (see ana_mutation_quality)


    %% Define MutQual and fixedmutation without outgroup isolates for SNP analysis


    [MutQual, ~] = ana_mutation_quality( Calls(:,~in_outgroup),Quals(:,~in_outgroup) );  % assumes quals already inverted
    fixedmutation = fixedmutation_tree; % fixed mutation initialized as the same as before, BUT...
    fixedmutation(MutQual < 1,:) = 0; % ...and not where there is only one type of non-N call


    %% remove maNT where maf is too low


    maNT( maf < params.min_maf_for_analysis ) = 0; % less stringent here


    %% Define goodpos as positions where there is at least one fixed mutation across all samples in the cluster


    goodpos =  sum( fixedmutation,2 ) > 0 ;
    fprintf(1,['Number of positions with fixed mutations: ' num2str(sum(goodpos)) '\n']);


    %% find recombinant SNPs. IF all positions are removed, continnue


    if sum(goodpos)>1
        involved_in_non_snp_event = findrecombinantSNPs(counts(:,:,~in_outgroup),anc_nti,p,goodpos);
        goodpos = goodpos&~involved_in_non_snp_event;
    end

    P = sum(goodpos);
    if P==0
        continue
    end

    %% get annotations


    annotations_cluster = annotate_mutations_gb_lineage_assembly(REFGENOMEFOLDER,p(goodpos), anc_nti(goodpos));
    annotations_full = append_annotations_jsb(annotations_cluster, anc_nti(goodpos), Calls(goodpos,:), counts(:,goodpos,:), params.promotersize);

    N = size(annotations_full,2);
    [annotations_full(1:N).cluster] = deal(CladeNumer);


    %% Add samplenames and mutation boolean to annotations


    has_mutation_ingroup_goodpos_no_reco = Calls(goodpos,~in_outgroup)>0&Calls(goodpos,~in_outgroup)~=repmat(anc_nti(goodpos),1,Nsamples_ingroup);

    MutationSampleNames = arrayfun(@(x) {SampleNames_ingroup(has_mutation_ingroup_goodpos_no_reco(x,:)==1)},(1:P)');
    MutationHasMutation= arrayfun(@(x) {has_mutation_ingroup_goodpos_no_reco(:,has_mutation_ingroup_goodpos_no_reco(x,:)==1)},1:P);
    for n = 1:N
        annotations_full(n).IsolateNamesWithMutation=MutationSampleNames(n);
        annotations_full(n).muts_hasmutation=MutationHasMutation(n);
    end

    annotations_full = fillemptyannotations(annotations_full);

    %% Get calls for analysis


    Calls_for_analysis=maNT(goodpos,:); % Note that this includes the outgroup
    Calls_for_analysis( maf(goodpos,:)      < params.min_maf_for_analysis ) = 0;
    Calls_for_analysis( coverage(goodpos,:) < params.min_coverage_for_analysis ) = 0;


    %% get tree directory


    cluster_tree_dir = [char(cladedir) '/' 'Tree'];

    % Directory for tree
    if ~exist(cluster_tree_dir,'dir')
        mkdir(cluster_tree_dir)
        copyfile('dnapars', [cluster_tree_dir '/dnapars'])
        copyfile('dnapars.app', [cluster_tree_dir '/dnapars.app'])
    end


    %% get calls for tree, true where different than anscestor

    tree_snps = Calls_for_analysis(:,~in_outgroup)~=0 & Calls_for_analysis(:,~in_outgroup) ~= repmat(anc_nti(goodpos),1,sum(~in_outgroup));
    tree_snps = Calls_for_analysis~=0 & Calls_for_analysis ~= repmat(anc_nti(goodpos),1,Nsample);

    %% Samplenames for tree and nucleotides for running DNApars

    % ancestral nucleotides
    calls_for_tree = [Calls_for_analysis(:,~in_outgroup) anc_nti(goodpos)];
    calls_ingroup_samples=zeros(size(calls_for_tree));
    calls_for_tree(calls_for_tree>0)=params.NTs(calls_for_tree(calls_for_tree>0));
    calls_for_tree(calls_for_tree==0)='N';
    names_for_tree = [SampleNames_ingroup {'anscestral_nucleotides'}];

    %% make phylogenies

    cluster_tree_dir = [char(cladedir) '/' 'Tree'];
    fprintf(1,'Making a tree...\n');
    % Directory for tree
    cluster_tree_dir = [char(cladedir) '/' 'Tree'];
    if ~exist(cluster_tree_dir,'dir')
        mkdir(cluster_tree_dir)
        copyfile('dnapars', [cluster_tree_dir '/dnapars'])
        copyfile('dnapars.app', [cluster_tree_dir '/dnapars.app'])
    end
    % make phylogeny
    cd(cluster_tree_dir);
    generate_parsimony_tree(calls_for_tree, names_for_tree, subjects_cluster_ingroup);
    cd(workingdir)


    %% add data to clades structure


    % clades.included_after_strict_filtering{k} = include;
    clades.anc_nti{k} = anc_nti;
    clades.goodpos{k}=goodpos;
    clades.outgroup{k}=in_outgroup;
    clades.calls_analysis{k}=Calls_for_analysis;
    clades.samplenames{k}=SampleNames;
    clades.p{k}=p;
    clades.has_mutation_ingroup_goodpos_no_reco{k}=has_mutation_ingroup_goodpos_no_reco;
    clades.annotations{k}=annotations_full;
    clades.tree_snps{k}=tree_snps;
    clades.cladenumber(k) = str2double(regexprep(CladeName,'.+clade_',''));

    
end
end