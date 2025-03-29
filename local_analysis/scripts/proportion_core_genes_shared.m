function proportion_core_genes_shared(uclusters,DM,str)
%% load roary structure


if str=="cacnes"
    PresenceAbsense = readtable('data/pangenome_cacnes_lineages/gene_presence_absence.csv');
    cladestr="cacnes_clade_";
    
elseif str=="sepi"
    PresenceAbsense = readtable('data/pangenome_sepi_lineages/gene_presence_absence.csv');
    cladestr="sepi_clade_";
end

%% parse columns 

Columns = string(PresenceAbsense.Properties.VariableNames);
keep=startsWith(Columns,cladestr);
C_table = Columns(keep);
PresenceAbsense=table2array(PresenceAbsense(:,keep));

%% re-index columns to match list of cluster numbers C
C_table=str2double(strrep(C_table,cladestr,''));

shared_both = intersect(uclusters,C_table);
keep=ismember(uclusters,shared_both);
DM=DM(keep,keep);
uclusters=uclusters(keep);

idx=arrayfun(@(x) find(x==C_table), uclusters);
PresenceAbsense=PresenceAbsense(:,idx);

%% turn PresenceAbsense into array which is true where there is a gene

Present=~cellfun(@isempty,PresenceAbsense);
PercentCladesWithGene=mean(Present,2);

percent_cutoffs = 0:.05:1;
P = numel(percent_cutoffs);
jaccard_d = cell(P,1);
number_same = cell(P,1);

for p =1:P
    meets_cutoff = PercentCladesWithGene>=percent_cutoffs(p);
    genes = Present(meets_cutoff,:);

    jaccard_d{p}=jaccard(genes');

    [~,C]=size(genes);
    number_shared = zeros(C);
    for i = 1:C
        for j = 1:C
            if i~=j
                number_shared(i,j)=sum(sum(genes(:,[i j]),2)==2);
            end
        end
    end
    number_same{p}=number_shared;
end

%% plot

for p =P:-1:1
    X = DM(:);
    Y_number = number_same{p}(:);
    scatter(X,Y_number,'filled','MarkerFaceColor',1.-([percent_cutoffs(p) percent_cutoffs(p) percent_cutoffs(p)]),'MarkerEdgeColor','black')
    hold on

end

n_core = sum(PercentCladesWithGene==1);

plot([0 max(X)],[n_core n_core],'LineWidth',1,'LineStyle','--')
ylim([1400 2400])
xlim([0 4.5*10^4])
