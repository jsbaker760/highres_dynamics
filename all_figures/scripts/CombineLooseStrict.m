function Combined = CombineLooseStrict(Strict,Loose,tbl)
%% initialize data


Combined = zeros(size(Strict));


%% iterate through subjects


uSID = unique(tbl.SID);

for u = 1:numel(uSID)
    % true where subject is subject
    bool=tbl.SID==uSID(u);
    if sum(bool)>0
        % slices for this subject
        ATight = Strict(bool,:);
        ALoose= Loose(bool,:);
        % true if ever found at any TP for strict dats
        FoundTightAny=repmat(sum(ATight>0,1)>0,size(ATight,1),1);
        % Missing
        MissingTight=ATight==0&FoundTightAny;
        %
        ATight(MissingTight)=ALoose(MissingTight);
        ATight(sum(ATight,2)>1,:)=ATight(sum(ATight,2)>1,:)./sum(ATight(sum(ATight,2)>1,:),2);
        Combined(bool,:)=ATight;
    end
end