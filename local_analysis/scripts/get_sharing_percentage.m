function get_sharing_percentage(M,LineagesPresent,str)

UniqueSubjects = unique(M.SID);
ever_found_sepi = arrayfun(@(x) {sum(LineagesPresent(x==M.SID,:),1)>0}, UniqueSubjects);
ever_found_sepi=vertcat(ever_found_sepi{:});

rm = sum(ever_found_sepi,2)==0;
UniqueSubjects(rm)=[];
ever_found_sepi(rm,:)=[];

Families = get_family_numbers(UniqueSubjects);
UniqueFamilies = unique(Families);
NFoundFamily = arrayfun(@(x) {sum(ever_found_sepi(x==Families,:),1)}, UniqueFamilies);
NFoundFamily=vertcat(NFoundFamily{:});

TotalFound = sum(sum(NFoundFamily,1)>0);
SharedInFamily=sum(sum(NFoundFamily>1,1)>0);

['Proportion ' str ' clades with any evidence of sharing from either data source amongst family members =' char(string(SharedInFamily/TotalFound))]