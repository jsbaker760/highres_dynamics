function D=get_dm_flexible_genecontent(C,str,cutoff_core)
%% load roary structure

if str=="cacnes"
    PresenceAbsense = readtable('data/Roary/cacnes_homologues_roary90/gene_presence_absence.csv');
    cladestr="cacnes_clade_";
    
elseif str=="sepi"
    PresenceAbsense = readtable('data/Roary/sepi_homologues_roary90/gene_presence_absence.csv');
    cladestr="sepi_clade_";
end

%% parse columns 

Columns = string(PresenceAbsense.Properties.VariableNames);
keep=startsWith(Columns,cladestr);
C_table = Columns(keep);
PresenceAbsense=table2array(PresenceAbsense(:,keep));

%% re-index columns to match list of cluster numbers C
C_table=str2double(strrep(C_table,cladestr,''));
idx=arrayfun(@(x) find(x==C_table), C);
PresenceAbsense=PresenceAbsense(:,idx);

%% turn PresenceAbsense into array which is true where there is a gene

Present=~cellfun(@isempty,PresenceAbsense);

PercentCladesWithGene=mean(Present,2);


Present(PercentCladesWithGene>=cutoff_core,:)=[];

%% make distance matrix from gene presence

D=jaccard(Present');
