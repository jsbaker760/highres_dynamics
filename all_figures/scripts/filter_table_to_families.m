function Mfilt = filter_table_to_families(M)
% 
AssignedReads=sum(M.BrackenReadsAssignedGenus,2);

nreadscutoff=50000;
EnoughCoverage = AssignedReads>nreadscutoff&~(sum(M.BrackenSpecies,2)==0)&~(sum(M.BrackenGenus,2)==0)
f=figure;hold on
histogram(AssignedReads,'BinWidth',10000)
xlim([0 1000000])
plot([nreadscutoff nreadscutoff],[0 40],'Color','red','LineStyle','--','LineWidth',1)
 MEnoughCoverage=M(EnoughCoverage,:)

 text(nreadscutoff,40,'cutoff=50000 reads')

saveas(f,'ManuscriptFigures/Supplemental_qc_tubes_removed_lowreads.svg')


[uSID,ai,ic]=unique(MEnoughCoverage.SID)

FamilyNumbers = get_family_numbers(uSID)

uFamilyNumbers = unique(FamilyNumbers);

HasBoth = arrayfun(@(x) any(contains(uSID(FamilyNumbers==x),"P"))&any(~contains(uSID(FamilyNumbers==x),"P")), uFamilyNumbers)
remove =ismember(get_family_numbers(MEnoughCoverage.SID),uFamilyNumbers(~HasBoth))|isnan(get_family_numbers(MEnoughCoverage.SID))|isnan(MEnoughCoverage.Age)'

tubes_removed1=(MEnoughCoverage(remove,:));
Mfilt = MEnoughCoverage(~remove,:);


F  = get_family_numbers(Mfilt.SID); 
U = unique(F);
Nperfam=arrayfun(@(x) numel(unique(Mfilt.SID(F==x))) , U);

rm2=ismember(F,U(Nperfam<2));
tubes_removed2=(Mfilt(rm2,:))

Mfilt(rm2,:)=[]


tubes_removed=[tubes_removed1 ; tubes_removed2];
writetable(tubes_removed(:,1:6),'ManuscriptFigures/Supplemental_qc_tubes_removed_no_family_members.csv')


end