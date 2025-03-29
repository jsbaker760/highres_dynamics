function is_consistent =  check_consistency(DM,clusters)

%% returns true if and only if all isolates are at least as far from intercluster isolates as intracluster isolates
is_consistent = true;
clustered =clusters>0;
DM=DM(clustered,clustered');
clusters = clusters(clustered);

C = max(clusters);
for c = 1:C
    isc = clusters==c;
    if any(max(DM(isc',isc),[],2) > min(DM(isc,~isc),[],2))
        is_consistent = false;
        break
    end
end