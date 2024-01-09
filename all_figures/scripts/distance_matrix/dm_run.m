function [distance_matrix] = dm_run(Calls_dist)
dbstop if error
Calls_dist=gpuArray(int8((Calls_dist)));
numsamples=gpuArray(int16(size(Calls_dist,2)));
goodcall = gpuArray(Calls_dist>0);
distance_matrix=gpuArray(single(zeros(numsamples)));
tic;
for s=1:numsamples
    sample_calls = repmat(Calls_dist(:,s),1,numsamples);
    distance_matrix(s,:)=sum(Calls_dist~=sample_calls & goodcall & sample_calls>0);
end
toc
distance_matrix=gather(distance_matrix);
distance_matrix(eye(size(distance_matrix))==1)=0;
end