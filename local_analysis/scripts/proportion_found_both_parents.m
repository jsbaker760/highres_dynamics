function proportion_found_both_parents(M,Abundances,str)

bool=sum(Abundances,2)>0;
M=M(bool,:);
Abundances(~bool,:)=[];

%% bring to individual parents

uSID = unique(M.SID);
uSID(~contains(uSID,"P"))=[];
Abundances = arrayfun(@(x) {sum(Abundances(x==M.SID,:),1)>0}, uSID);
Abundances=vertcat(Abundances{:});


%% remove one family things

FamilyNumber=get_family_numbers(uSID);
n_in_family=arrayfun(@(x) sum(FamilyNumber==x), FamilyNumber);
Abundances(n_in_family<2,:)=[];
FamilyNumber(n_in_family<2)=[];

%% remove non-two parent things


UniqueFamilyNumbers = unique(FamilyNumber);
NFoundFamily = arrayfun(@(x) {sum(Abundances(x==FamilyNumber,:),1)}, UniqueFamilyNumbers);
NFoundFamily=vertcat(NFoundFamily{:});

TotalFound = sum(sum(NFoundFamily,1)>0);
SharedInFamily=sum(sum(NFoundFamily>1,1)>0);

NotSharedInFamily = TotalFound-SharedInFamily;

[char(string(NotSharedInFamily)) '/' char(string(TotalFound)) str ' lineages are unshared']