function [MinimumNumberTransmissions,Transmissions,AllPossibleTransmissions]=get_all_transmission_counts_pairwise(L,min_subject,min_other)
%% Make a table of all possible transmissions

Transmissions=table;
idx = 1;
for i =1:size(L,1)
    Subjects=char(L.samplenames{i}(~L.outgroup{i}));
    Subjects=string(Subjects(:,1:3));
    TreeSNPs=L.tree_snps{i}';
    uSubjects = unique(Subjects);
    for s = 1:numel(uSubjects)
        source_subject=uSubjects(s);
        other_subjects=uSubjects(~ismember(uSubjects,source_subject));
        numtransmissions=zeros(numel(other_subjects),1);
        numisolates=zeros(numel(other_subjects),1);
        compared=false(size(numtransmissions));
        for o = 1:numel(other_subjects)
            recipient = other_subjects(o);
            is_source=Subjects==source_subject;
            is_recipient=Subjects==recipient;
            Nsource=sum(is_source);
            Nrecipient = sum(is_recipient);
            if Nsource>min_other&Nrecipient>min_subject
                rows_include=is_recipient|is_source;
                Tree = TreeSNPs(rows_include,:);
                Tree(:,sum(Tree,1)==size(Tree,1))=[];
                is_recipient=is_recipient(rows_include);
                numtransmissions(o)=count_transmission_on_trees(Tree,is_recipient);
                compared(o)=true;
                numisolates(o)=sum(is_recipient);
            end
        end
        if any(compared)
            Transmissions.Source(idx)=source_subject;
            Transmissions.Recipients{idx}=other_subjects(compared);
            Transmissions.Ntransmissions{idx}=numtransmissions(compared);
            Transmissions.LineageNumber(idx)=L.cladenumber(i);
            Transmissions.NisolatesRecipient{idx}=numisolates(compared);
            Transmissions.NisolatesSource{idx}=sum(ismember(Subjects,source_subject))';

            idx = idx+1;
        end

    end
end

%% plot the minimum
AllPossibleTransmissions=Transmissions;
Transmissions.TotalTransmissions=cellfun(@sum,Transmissions.Ntransmissions);
idxMinNperLineage=arrayfun(@(x) find(Transmissions.LineageNumber==x&Transmissions.TotalTransmissions==min(Transmissions.TotalTransmissions(Transmissions.LineageNumber==x)),1), unique(Transmissions.LineageNumber));
Transmissions=Transmissions(idxMinNperLineage,:);

%
MinimumNumberTransmissions=vertcat(Transmissions.Ntransmissions{:});
numisolates=vertcat(Transmissions.NisolatesRecipient{:});
MinimumNumberTransmissionsLineage = arrayfun(@(x) sum(~L.outgroup{L.cladenumber==x}), Transmissions.LineageNumber);

histogram(MinimumNumberTransmissions,'BinWidth',1,'BinEdges',0:1:12);
xlabel('N genotypes transmitted (most parsimonious)')
ylabel('N instances')
xticks(1.5:1:12.5)
xticklabels(1:1:12)
xlim([1 12])
pbaspect([1 1 1])

%%
