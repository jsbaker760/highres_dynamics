function coefs = ACcoef(DM,E)

%% calculate local clustring coefficients
G = DM<=E;
G(eye(size(G))==1)=0;

n=length(G);
coefs=(single(zeros(n,1)));


kall =sum(G);
kall = (kall.^2)-kall;
for u=1:n
    k = kall(u);
    if k 
        V=G(u,:);
        V=G(V,V);
        coefs(u)=sum(V(:))/k;
    end
    
end
