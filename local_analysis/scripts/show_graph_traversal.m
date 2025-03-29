function show_graph_traversal(DMin,clustersin)

DM=DMin;
clusters=clustersin;
W=20;

figure
A=DM<=W;
G= graph(A,'omitselfloops');
H = plot(G,'Layout','subspace');
title('connectivity of all isolates, including those which were excluded')
% add colors

edge_clusters = clusters(G.Edges.EndNodes);
NodeI = edge_clusters(:,1);
NodeJ = edge_clusters(:,2);
intracluster = NodeI==NodeJ&NodeI>0&NodeJ>0;
cluster_to_bad = NodeI~=NodeJ&(NodeI>0|NodeJ>0)&(NodeI<=0|NodeJ<=0);
intercluster= (NodeI~=NodeJ)&(NodeI>0&NodeJ>0);
% singleton_to_singleton = NodeI<0&NodeJ<0;

E = numel(cluster_to_bad);
N = size(A,1);

EdgeCData = zeros(E,3);
NodeCData = zeros(N,3);

NodeCData(clusters<=0,1)=1;

EdgeCData(intracluster,:)=0;
EdgeCData(cluster_to_bad,1)=1;
EdgeCData(intercluster,3)=1;
% EdgeCData(singleton_to_singleton,2)=1;

H.EdgeColor=EdgeCData;
H.NodeColor=NodeCData;

H.LineWidth=.5;

%% example: everythign traversible by 1012

bins = conncomp(G,'OutputForm','cell');
c = bfsearch(G,1012);

allnodes =  [1012; c];
DM = DM(allnodes,allnodes);
clusters=clusters(allnodes);
%%
figure
A=DM<=W;
G= graph(A,'omitselfloops');
names = arrayfun(@(x) strjoin(["cluster" string(clusters(x)) "node" string(x)],' '), 1:size(A,1))

H = plot(G,'Layout','force','NodeLabel',names);
%H = plot(G,'Layout','force','UseGravity',true);
% add colors

edge_clusters = clusters(G.Edges.EndNodes);
NodeI = edge_clusters(:,1);
NodeJ = edge_clusters(:,2);
intracluster = NodeI==NodeJ&NodeI>0&NodeJ>0;
cluster_to_bad = (NodeI<=0|NodeJ<=0);
intercluster= (NodeI~=NodeJ)&(NodeI>0&NodeJ>0);
% singleton_to_singleton = NodeI<0&NodeJ<0;

E = numel(cluster_to_bad);
N = size(A,1);

EdgeCData = zeros(E,3);
NodeCData = zeros(N,3);

NodeCData(clusters<=0,1)=1;

EdgeCData(intracluster,:)=0;
EdgeCData(cluster_to_bad,1)=1;
EdgeCData(intercluster,3)=1;
% EdgeCData(singleton_to_singleton,2)=1;

H.EdgeColor=EdgeCData;
H.NodeColor=NodeCData;

H.LineWidth=.5

title('C. acnes isolates which were not clustered are connected to multiple self consistent lineage clusters')

%%
 

%%

clustered=clusters>0

DM = DM(clustered,clustered);
clusters=clusters(clustered);

figure
A=DM<=W;
G= graph(A,'omitselfloops');
names = arrayfun(@(x) strjoin(["cluster" string(clusters(x)) "node" string(x)],' '), 1:size(A,1))

%H = plot(G,'Layout','subspace','NodeLabel',names);
H = plot(G,'Layout','subspace','NodeLabel',names);
% add colors
title('exclusion of isolates allows self-consistency')

edge_clusters = clusters(G.Edges.EndNodes);
NodeI = edge_clusters(:,1);
NodeJ = edge_clusters(:,2);
intracluster = NodeI==NodeJ&NodeI>0&NodeJ>0;
cluster_to_bad = NodeI~=NodeJ&(NodeI>0|NodeJ>0)&(NodeI<=0|NodeJ<=0);
intercluster= (NodeI~=NodeJ)&(NodeI>0&NodeJ>0);
singleton_to_singleton = NodeI<0&NodeJ<0;

E = numel(cluster_to_bad);
N = size(A,1);

EdgeCData = zeros(E,3);
NodeCData = zeros(N,3);

NodeCData(clusters<=0,1)=1;

EdgeCData(intracluster,:)=0;
EdgeCData(cluster_to_bad,1)=1;
EdgeCData(intercluster,3)=1;
EdgeCData(singleton_to_singleton,2)=1;

H.EdgeColor=EdgeCData;
H.NodeColor=NodeCData;

H.LineWidth=1