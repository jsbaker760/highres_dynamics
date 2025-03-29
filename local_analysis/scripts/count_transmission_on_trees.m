function numtransmissions = count_transmission_on_trees(genotypes, is_recipient)

%Each row in genotypes is a genotype, each column is a mutation
%Each row in compartmenttimes is a genotype, each column is a compartment,
%values are the number of times it was found in that compartment

%Tami Lieberman, July 2015






%

genotypes_recipient=genotypes(is_recipient,:);

% found_on_recipient = arrayfun(@(x) sum(sum(sum(genotypes(x,:)==genotypes_recipient)==0,2)==0), 1:size(genotypes,1))'
genotypes_in_destination=is_recipient

% intialize
numtransmissions=0;

mutations_in_source=find(sum(genotypes(~is_recipient,:),1)>0);
    
% genotypes_in_destination=compartmenttimes(:,compartment_number_of_destination)>=mintimes;
mutations_in_recipient=find(sum(genotypes_recipient,1)>0);

unaccounted_for_shared_mutations=mutations_in_recipient(ismember(mutations_in_recipient,mutations_in_source));

%
% 
genotypes_containing_shared_mutation=sum(genotypes(:,unaccounted_for_shared_mutations),2)>0;
% 
% %% deal with special cases 
% % 
ismultiplestraininfection=false;
if ~ismultiplestraininfection
    %(1) ancestor of patient found in organ
    if sum(sum(genotypes(genotypes_in_destination,:),2)==0)>0 %ancestor in organ
        numtransmissions=numtransmissions+1;
   %(2) where the ancestor isn't found in organ but a unique descendent of it is
    elseif  sum(genotypes_in_destination & ~genotypes_containing_shared_mutation) > 0
        numtransmissions=numtransmissions+1;
    %(3) no genotype, including ancestor found more than
    %MIN_OBSERVATIONS_FOR_TRANSMISSION_ANALYSIS -- infer a single
    %transmission event of ancestor
    % (alternatively this could have been chosen to be an excluded case)
    elseif sum(genotypes_in_destination) == 0 
        numtransmissions=numtransmissions+1;
    end
end


% iteratively account for mutations

% SN={}
while numel(unaccounted_for_shared_mutations)>0

    numtransmissions=numtransmissions+1;
    %find genotypes that contribute to these
    %unshared mutations
    %pick genotypes with fewest mutations first
    %first look for genotypes in organ 
    
     genotypes_containing_shared_mutation=find(sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_destination));
%         SN{end+1}=SampleNames(genotypes_containing_shared_mutation);
    if isempty(genotypes_containing_shared_mutation)
        error(1,'did not find a genotype with unaccounted for mutation in destination organ');
    end
    
    minmuts=min(sum(genotypes(genotypes_containing_shared_mutation,:),2));
    mingenotype=genotypes_containing_shared_mutation(find(sum(genotypes(genotypes_containing_shared_mutation,:),2)==minmuts,1));

    %remove all unaccounted for mutations carried by this minimal genotype
    mutations_remove=ismember(unaccounted_for_shared_mutations,find(genotypes(mingenotype,:)));
%     SampleNames(mingenotype)
    unaccounted_for_shared_mutations(mutations_remove)=[];
    
    
end


%%
% SN=[ {SampleNames(genotypes_in_destination)} SN]
% SN = arrayfun(@(x) SN{x}(~ismember(SN{x},SN{x+1})), 1:(numel(SN)-1), 'UniformOutput', false)


%% sanity check

if numtransmissions==0
    fprintf(1,'Warning: no transmissions');
end