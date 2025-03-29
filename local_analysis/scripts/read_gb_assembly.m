function CDS = read_gb_assembly(RefGenomeDirectory)


% May need revisions for your reference genome

%%
%Edited by TDL 2018_08_13 to not save genenum -- can be generated quickly
% with genomic_position_all
% FMK: CDS generation using genbankFile.Sequence if exists rather than full
% genome. Previously led to wrong seq/translation for chr2 and onwards.
% 2019/11 JSB and TDL changed extractFeature sub-function of
% featurelocation div_add_in_nonCDS

nts='atcg';
rc='tagc';

%% Intialize
S=strsplit(RefGenomeDirectory,'/');
cladename = S(end);
% Load fasta
fastaFilename = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFilename);

% load genbank file
gbfilename = [RefGenomeDirectory '/' cladename{:} '.gbff'];
    if ~isfile(gbfilename)
            gbfilename = [RefGenomeDirectory '/prokka_out.faa'];
    end

    genbankFile = genbankread(gbfilename);


[ChrStarts, glength, ~, ScafNames]= genomestats(RefGenomeDirectory);

%%
tic; fprintf(1,'Reading in gb file...\n ');

L = numel(genbankFile);
S = numel(ScafNames);
clear CDS
CDS=cell(S,1);



parfor i=1:L
    
    ScafNames_i=ScafNames{i};
    
    fprintf(1,[ScafNames_i '...\n']);
    
    % FIND ANNOTATION SOURCE FOR SCAFFOLD SEQUENCE (GB or FASTA)
    
    % LOAD GENBANK
 
    scafSeq = genbankFile(i).Sequence;
    
    
    if isfield(genbankFile,'CDS') || isfield(genbankFile,'Features')
        genes = locustag_from_text(genbankFile(i).CDS) ;
        genes = div_add_in_nonCDS(genes, genbankFile(i).Features); %fixes many things that matlab's genbank reader fails to import
        if isfield(genbankFile, 'Sequence')
            CDS{i} = parse_all_locations_gb(genes, scafSeq); %#ok<*PFOUS>
        end
            %sort by position on genome
        if ~isempty(CDS{i})
            [~,sortedpositions]=sort([CDS{i}.loc1]);
            CDS{i}=CDS{i}(sortedpositions);
        end
    end
    
    
end


save([RefGenomeDirectory '/cds_sorted'],'CDS')

return
end