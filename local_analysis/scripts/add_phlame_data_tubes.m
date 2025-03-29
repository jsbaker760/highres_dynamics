function Mtubes=add_phlame_data_tubes(Mtubes,M)

%% get the indices for each unique tube

R = size(M,1);
IndividualTubes = arrayfun(@(x) strjoin([M.SID(x) string(M.TP(x)) M.location(x)]) , 1:R);
[uIndividualTubes,ia,ic]=unique(IndividualTubes);
index = arrayfun(@(x) {find(IndividualTubes==x)}, uIndividualTubes)';

%% For each type of abundance, get the one with the best assignment 
Mtubes.CacnesLineageAbundance=add_best_instance(M.CacnesLineageAbundance,index);
Mtubes.CacnesPhylotypeAbundance=add_best_instance(M.CacnesPhylotypeAbundance,index);
Mtubes.SepiLineageAbundance=add_best_instance(M.SepiLineageAbundance,index);
Mtubes.SepiPhylotypeAbundance=add_best_instance(M.SepiPhylotypeAbundance,index);
Mtubes.SepiPhylotypeAbundance=add_best_instance(M.SepiPhylotypeAbundance,index);



function best_instance=add_best_instance(Abundances, index)

best_instance = cell(numel(index),1);

for i = 1:numel(best_instance)
    instances = Abundances(index{i},:);
    assigned = sum(instances,2);
    best = max(assigned);
    if best==0
        best_instance{i}=zeros(1,size(instances,2));
    else
        best = instances(assigned==best,:);
        if size(best,1)>1 % of there's multiple instances with complete assignment
            best=best(1,:);
        end
        best_instance{i}=best;
    end
end

best_instance=vertcat(best_instance{:});


