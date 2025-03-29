function [f, f2]=unshared_clades_by_group(MsubjectTime,MinAssigned,MinAbundance,Abundances,ManuscriptColors)
%% remove samples without enough assignment and also low abundance
MsubjectTime_filt=MsubjectTime(sum(Abundances,2)>MinAssigned,:);
AbundancesFilt=Abundances(sum(Abundances,2)>MinAssigned,:);
AbundancesFilt(AbundancesFilt<MinAbundance)=0;

%% get family number of each subject and remove rows where there aren't at least two different subjects with enough assignment
FamilyNumbers = get_family_numbers(MsubjectTime_filt.SID);
bool=arrayfun(@(x) numel(unique(MsubjectTime_filt.SID(FamilyNumbers==x)))>1, FamilyNumbers);
MsubjectTime_filt=MsubjectTime_filt(bool,:);
AbundancesFilt=AbundancesFilt(bool,:);
FamilyNumbers=FamilyNumbers(bool);

%%  get all types ever found on subject
[UniqueSubjects, iu,ic]= unique(MsubjectTime_filt.SID);

IsEverFoundOnSubject = arrayfun(@(x) {sum(AbundancesFilt(x==MsubjectTime_filt.SID,:),1)>0}, UniqueSubjects);
IsEverFoundOnSubject=vertcat(IsEverFoundOnSubject{:});
%% get boolean which is true where this clade is never found on another family member

IsUniqueToSubjectInFamily=arrayfun(@(x) {IsEverFoundOnSubject(x,:)>0&sum(IsEverFoundOnSubject((1:numel(UniqueSubjects)~=x)&(FamilyNumbers(iu)==FamilyNumbers(iu(x))),:),1)==0} ,1:numel(UniqueSubjects));
IsUniqueToSubjectInFamily=vertcat(IsUniqueToSubjectInFamily{:});


%% get basic data structures for plotting

unshared_clade_abundances = AbundancesFilt(IsUniqueToSubjectInFamily(ic,:)&AbundancesFilt>0);
shared_clade_abundances = AbundancesFilt(~IsUniqueToSubjectInFamily(ic,:)&AbundancesFilt>0);
nclades_total = sum(IsEverFoundOnSubject,2);
nclades_unshared = sum(IsUniqueToSubjectInFamily,2);

%% plot the relative abundances of things found on other family members or not

f=figure;

X = MsubjectTime_filt.subject_cutotype(iu);
Y1=sum(nclades_unshared,2);
Y2=sum(nclades_total,2);
Y3=Y1./Y2;
subplot(2,1,1)
swarmchart(X,Y1,'filled','CData',ManuscriptColors.CutotypeColors(X,:))
xlim([0 4])
xticks(1:3)
xticklabels(["FC1-Children" "FC2-Children" "Parents"])
ylabel('Total number of unshared clades across timepoints');
yticks([1:10])
pv = ranksum(Y1(X==1),Y1(X==3));
text(1.5,.9,['p=' char(string(pv))])
pbaspect([1 1 1])

X = [ones(size(unshared_clade_abundances)) ; 2*ones(size(shared_clade_abundances))];
Y = [unshared_clade_abundances ; shared_clade_abundances];

subplot(2,1,2)

swarmchart(X,Y,'filled','o','MarkerFaceColor','black','MarkerEdgeColor','black')
hold on
xticks(1:2)
xticklabels(["Lineage never found on related subject","Lineage ever found on related subject"])
ylim([0 1])

pv = ranksum(Y(X==1),Y(X==2));
text(1.5,.9,['p=' char(string(pv))])
pbaspect([1 1 1])
end