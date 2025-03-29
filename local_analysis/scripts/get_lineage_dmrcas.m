function LineageDMRCAs= get_lineage_dmrcas(LineageData,M)


LineageDMRCAs = table;
idx = 1;
for i =1:size(LineageData,1)
    Subjects = LineageData.samplenames{i}(~LineageData.outgroup{i});
    Subjects=char(Subjects);
    Subjects=string(Subjects(:,1:3));
    SID = unique(Subjects);
    if any(startsWith(SID,"0"))
        error('shoudnt be including control subjects')
    else
        N = numel(SID);
        if isempty(LineageData.p{i})

        end
        TreeSNPs=LineageData.tree_snps{i};
        if isempty(TreeSNPs)
            TreeSNPs=zeros(1,numel(Subjects));
        end
        for n=1:N
            % basic information
            LineageDMRCAs.SID(idx)=SID(n);
            LineageDMRCAs.Species(idx)=LineageData.SpeciesName(i);
            LineageDMRCAs.LineageNumber(idx)=LineageData.cladenumber(i);
            % true where isolate belongs to subject n of N
            bool = Subjects==SID(n);
            LineageDMRCAs.Nisolates(idx)=sum(bool);
            % root to tip distance of each isolate in clade, and then of
            % each unique genotype
            LineageDMRCAs.RTT_clade(idx)= {sum(TreeSNPs,1)};
            LineageDMRCAs.RTT_clade_ugenotypes(idx)= {sum(unique(TreeSNPs','rows'),2)};
            LineageDMRCAs.TreeSNPs{idx}=TreeSNPs;
            LineageDMRCAs.Subjects{idx}=Subjects;
            %
            subject_muts = TreeSNPs(:,bool);
            others_muts = TreeSNPs(:,~bool);
            unique_subject_muts = unique(subject_muts','rows')';
            %%
            LineageDMRCAs.UniqueSNPs_genotype{idx} = arrayfun(@(x) sum(unique_subject_muts(:,x)==1&sum(unique_subject_muts(:,1:size(unique_subject_muts,2)~=x),2)==0&sum(others_muts,2)==0), 1:size(unique_subject_muts,2));
            %
            LineageDMRCAs.dDMRCA(idx) = sum(sum(subject_muts,2)==size(subject_muts,2));
            LineageDMRCAs.MRCA_equals_root(idx)=~any(sum(subject_muts,2)==size(subject_muts,2))|all(subject_muts(:)==0);
            %
            LineageDMRCAs.UniqueSNPs_subject_isolates{idx} = sum(subject_muts(sum(subject_muts,2)>0&sum(others_muts,2)==0,:));
            LineageDMRCAs.UniqueSNPs_subject_ugenotypes{idx} = sum(unique_subject_muts(sum(unique_subject_muts,2)>0&sum(others_muts,2)==0,:));
            %
            LineageDMRCAs.RTT_subject_isolates{idx} = sum(subject_muts,1);
            LineageDMRCAs.RTT_subject_ugenotypes{idx} = sum(unique_subject_muts,1);
            %
            subject_muts(sum(subject_muts,2)==size(subject_muts,2),:)=[];
            unique_subject_muts(sum(unique_subject_muts,2)==size(unique_subject_muts,2),:)=[];
            %
            LineageDMRCAs.SubjectTreeSNPs_unique_genotypes{idx}=unique_subject_muts;
            LineageDMRCAs.SubjectTreeSNPs_isolates{idx}=subject_muts;
            %
            LineageDMRCAs.dMRCA_subject_isolates{idx} = sum(subject_muts,1);
            LineageDMRCAs.dMRCA_subject_ugenotypes{idx} = sum(unique_subject_muts,1);
            %
            if N==1
                LineageDMRCAs.Unshared(idx)=true;
            else
                LineageDMRCAs.Unshared(idx)=false;
            end
            idx = idx+1;
        end
    end
end
CutotypesSubject = M.subject_cutotype;
LineageDMRCAs.Cutotype = arrayfun(@(x) CutotypesSubject(M.SID==x), LineageDMRCAs.SID);
LineageDMRCAs.Age = arrayfun(@(x) M.Age(M.SID==x), LineageDMRCAs.SID);
end