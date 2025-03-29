function [PanarooDataCacnes, PanarooDataSepi]= parse_homolog_clusters(PanarooCacnes,PanarooSepi)

% get part of Panaroo output that has locustags for each cluster
[GeneArrayCacnes,ClusterNumbersCacnes,~]= parse_panaroo(readtable(PanarooCacnes));
[GeneArraySepi,ClusterNumbersSepi,~]= parse_panaroo(readtable(PanarooSepi));

% resize it to be table
PanarooDataCacnes=expand_Panaroo(GeneArrayCacnes,ClusterNumbersCacnes);
PanarooDataSepi=expand_Panaroo(GeneArraySepi,ClusterNumbersSepi);

end
%% resize to be 1D
function T=expand_Panaroo(GeneArray,CladeNumber)

[Nclusters,Nclades]=size(GeneArray);
% give a number to each row (one for each homologue)
% give a number to each column (repeats down for each lineage)
HomologueNumber=repmat((1:Nclusters)',1,Nclades);
CladeNumber=repmat(CladeNumber,Nclusters,1);

% resizing
GeneArray=GeneArray(:);
HomologueNumber=HomologueNumber(:);
CladeNumber=CladeNumber(:);

% remove bad rows
remove=cellfun(@isempty,GeneArray);
GeneArray(remove)=[];
HomologueNumber(remove)=[];
CladeNumber(remove)=[];
% number of rows
R = numel(GeneArray);
% split GeneArray strings where there's more than one;
GeneArray=arrayfun(@(x) {strsplit(x{:},';')},GeneArray);

% make cluster numbers a cell
Ngenes = cellfun(@numel,GeneArray);

% add extra numbers wheres there's more than one gene
HomologueNumber=arrayfun(@(x) {HomologueNumber(x)*ones(Ngenes(x),1)}, 1:R);
CladeNumber=arrayfun(@(x) {CladeNumber(x)*ones(Ngenes(x),1)}, 1:R);
% expand

HomologueNumber=vertcat(HomologueNumber{:});
GeneArray=horzcat(GeneArray{:})';
CladeNumber=vertcat(CladeNumber{:});

T=table;
T.HomologueNumber=HomologueNumber;
T.GeneTag=GeneArray;
T.CladeNumber=CladeNumber;      

end

%%
function [GeneArray,ClusterNumbers,HomologueInfo]=parse_panaroo(PanarooOutput)

Columns = string(PanarooOutput.Properties.VariableNames);
keep=contains(Columns,"_clade_");

GeneArray=table2array(PanarooOutput(:,keep));
HomologueInfo=PanarooOutput(:,~keep);

ClusterNumbers=Columns(keep);
ClusterNumbers=arrayfun(@(x) {strsplit(x,"_")}, ClusterNumbers);
ClusterNumbers=arrayfun(@(x) str2double(x{:}(3)), ClusterNumbers);
end
