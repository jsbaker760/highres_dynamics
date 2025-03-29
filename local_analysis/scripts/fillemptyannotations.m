function annotations_full = fillemptyannotations(annotations_full)


    allnames = {'gene_num','scaffold','pos','gene','protein','protein_id','strand','loc1','loc2','note','locustag','text','oldlocustag','Sequence','translation','nt_pos','aa_pos','codons','AA','NonSyn','gene2','locustag2','distance2','protein2','gene1','locustag1','distance1','protein1','WARNING','anc','AApos','annotation','type','AAs','nts','IsolateNamesWithMutation','cluster','muts_hasmutation','muts','locustag_CDHITcluster','locustag1_CDHITcluster','locustag2_CDHITcluster'};
    lacking = ~ismember(allnames, fieldnames(annotations_full));
    if any(lacking)
        doesnt_have = allnames(lacking);
        idx = 1:numel(annotations_full);
        for f = 1:numel(doesnt_have)
            annotations_full = arrayfun(@(idx) setfield(annotations_full(idx),char(doesnt_have(f)),[]),idx);
        end
    end
    % set gene_numt o 0 if it is empty
    for i = 1:numel(annotations_full)
    if isempty(annotations_full(i).gene_num)
        annotations_full(i).gene_num=0;
    end
    end
end