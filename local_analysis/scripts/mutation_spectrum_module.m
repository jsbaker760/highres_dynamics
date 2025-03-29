function [ mutationmatrix, mut_observed, typecounts, prob_nonsyn ] = ...
    mutation_spectrum_module( annotation_full, NTs)

% Initialize vectors
mutationmatrix=zeros(4,4); % ATCG -> ATCG matrix; gathered into mut_observed (6 entry vector)
typecounts=zeros(4,1); %NSPI % used to calculate prob_nonsyn

% Count mutations
for i=1:numel(annotation_full) % loop through positions at which there is a SNP
    anc = annotation_full(i).anc; % ancestor NT at this position
    new = annotation_full(i).nts; %all NTs found at this position
    if anc > 0 && sum(anc == new) % if ancestor is known and ancestor NT is in annotation_full
        anc = find(NTs==anc); % convert to numbers
        new = find(ismember(NTs,new)); % convert to numbers
        % Remove ancestor from new
        new = new( new ~= anc );
        % Remove N from new
        new = new( new ~= 0 );
        if numel(new)==0 % if nothing left
            fprintf(1,['Warning: No mutation found in Cluster ' num2str(annotation_full(i).cluster) ' position ' num2str(annotation_full(i).pos) '\n'])
        elseif numel(new)==1 % if one mutation
            % Update mutation matrix
            mutationmatrix(anc,new)=mutationmatrix(anc,new)+1;
            % Count type (use annotation_full)
            if annotation_full(i).type=='N' % nonsyn mut
                typecounts(1)=typecounts(1)+1;
            elseif annotation_full(i).type=='S' % syn mut
                typecounts(2)=typecounts(2)+1;
            elseif annotation_full(i).type=='P' % promoter mut
                typecounts(3)=typecounts(3)+1;
            elseif annotation_full(i).type=='I' % intergenic mut
                typecounts(4)=typecounts(4)+1;
            else
                fprintf(1,['Warning: Unrecognized mutations type ' annotation_full(i).type ' found in Cluster ' num2str(annotation_full(i).cluster) ' position ' num2str(annotation_full(i).pos) '.\n'])
            end
        elseif numel(new) > 1 % if more than one mutation
            fprintf(1,['Warning: Multiple mutations found in Cluster ' num2str(annotation_full(i).cluster) ' position ' num2str(annotation_full(i).pos) '.\n'])
            % Update mutation matrix
            for m=1:numel(new) % once for each mutation
                mutationmatrix(anc,new(m))=mutationmatrix(anc,new(m))+1;
            end
        end
    else
        fprintf(1,['Warning: Ancestor not found in Cluster ' num2str(annotation_full(i).cluster) ' position ' num2str(annotation_full(i).pos) '.\n'])
    end
end

% Count six types of mutations
mut_observed = div_matrix2_6types(mutationmatrix); % tally of mutations, NOT normalised

% Calculate fraction of nonsynonymous mutations (only N and S)
prob_nonsyn = typecounts(1)/(typecounts(1)+typecounts(2));


end
