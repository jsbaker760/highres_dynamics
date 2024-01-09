function [calls_for_tree_filtered, SampleNames_tagoutgroup_filtered,include]=filter_samples_for_strict_tree(calls_for_tree,SampleNames_tagoutgroup)
%% get samplenames and subject ID of each samplename


SN = SampleNames_tagoutgroup(~startsWith(SampleNames_tagoutgroup,"OUTGROUP"));
Subjects = arrayfun(@(x) string(regexp(x,'[0-8][APC][ABC]','match')) , SN);


%% get plate and well 


% get plate used for picking isolates
[pick_plate,~] = samplenames2pickingwells(SN);


%% for each subject, find out if there is enough evidence to be strongly confident that there really is subject's isolates in clade


uSubjects =unique(Subjects);
exclude = false(size(uSubjects));
for i = 1:numel(uSubjects)
    % get plates this subjects isolates on, plates unique to subject, n
    % isolates per subject
    plates_this_subject = unique(pick_plate(uSubjects(i)==Subjects));
    plate_unique_to_subject = zeros(size(plates_this_subject));
    n_isolates_per_subject  = zeros(size(plates_this_subject));
    for u = 1:numel(plate_unique_to_subject)
        n_isolates_per_subject(u) = sum(plates_this_subject(u)==pick_plate&uSubjects(i)==Subjects');
        if sum(plates_this_subject(u)==pick_plate&uSubjects(i)~=Subjects)==0
            plate_unique_to_subject(u)=true;
        end
    end
    % need to have at least 2 isolates and one on a plate nobody else in
    % the lineage has
    n_isolates_unique_plates = sum(n_isolates_per_subject(plate_unique_to_subject==1));
    exclude(i)= sum(n_isolates_per_subject<=2)&(sum(n_isolates_unique_plates)==0);
end

% remove subjects without extreme confidence
include = ~exclude;
include = find(ismember(Subjects,uSubjects(include)));
include = ismember(1:numel(SampleNames_tagoutgroup),include)|startsWith(SampleNames_tagoutgroup,"OUTGROUP");

% remove these subjects' names and positions which are no longer relavent
SampleNames_tagoutgroup_filtered=SampleNames_tagoutgroup(include);
calls_for_tree_filtered = calls_for_tree(:,include);
rmpos = arrayfun(@(x) numel(unique(calls_for_tree_filtered(x,:)))==1, 1:size(calls_for_tree_filtered));

calls_for_tree_filtered(rmpos,:)=[];