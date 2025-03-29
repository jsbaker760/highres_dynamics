function plot_sweeps_individual(SubjectClusterPairs,anc_nts,MRCAs,MRCA2)

%% initialize
C=numel(SubjectClusterPairs);
%% Iterate throught subject cluster pairs and find differences in DMRCA
interval_differences = cell(C,1);
for c = 1:C
    if ~isempty(MRCAs{c})&~isempty(MRCA2{c})
        MRCAs_1 = (MRCAs{c});
        MRCAs_2 = (MRCA2{c});
        anc_nt = anc_nts{c};
        SubjectClusterPairs(c)
        mutations_found_previous_timepoint =mean(MRCAs_1,2)~=0&mean(MRCAs_1,2)~=1;
        mutations_MRCA_last_timepoint = mean(MRCAs_2,2)==1;
        new = sum(mutations_found_previous_timepoint==0&mutations_MRCA_last_timepoint==1);
        interval_differences{c}=sum(new);
    end
end

%% pull out differences

idx = arrayfun(@(x) {repmat(x,numel(interval_differences{x}),1)} , 1:C);
idx = vertcat(idx{:});
interval_differences=horzcat(interval_differences{:});
instances = SubjectClusterPairs(idx);
has_difference = interval_differences>0;
nonzero_differences=interval_differences(has_difference);
instances_with_differences = instances(has_difference);

[~,i]=sort(nonzero_differences,'descend');

instances_with_differences=instances_with_differences(i);
nonzero_differences=nonzero_differences(i);

%% plot instances with difference

figure
subplot(2,1,1)
histogram(interval_differences(contains(instances,"cacnes")),'BinWidth',.9,'FaceColor','white')
xlabel('Number of timepoint specific SNPs not found at a previous timepoint')
ylabel('number of instances')
title('C. acnes')
xlim([0 35])
ylim([0 20])
subplot(2,1,2)
histogram(interval_differences(contains(instances,"sepi")),'BinWidth',.9,'FaceColor','black')
xlabel('Number of timepoint specific SNPs not found at a previous timepoint')
ylabel('number of instances')
title('S. epidermidis')
xlim([0 55])
ylim([0 20])