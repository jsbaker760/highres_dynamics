function D=get_dm_genecontent(C,PresenceAbsense,cladestr)
%% parse columns

Columns = string(PresenceAbsense.Properties.VariableNames);
keep=startsWith(Columns,cladestr);
C_table = Columns(keep);
PresenceAbsense=table2array(PresenceAbsense(:,keep));

%% re-index columns to match list of cluster numbers C
C_table=str2double(strrep(C_table,cladestr,''));
idx=arrayfun(@(x) find(x==C_table), C);
PresenceAbsense=PresenceAbsense(:,idx);

%% turn PresenceAbsense into array which is true where there is a gene
Present=~cellfun(@isempty,PresenceAbsense);
% remove empties
Present(sum(Present,2)==0,:)=[];


%% make distance matrix from gene presence
D=jaccard(Present');

end
