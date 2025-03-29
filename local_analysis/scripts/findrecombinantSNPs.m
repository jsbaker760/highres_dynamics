function involved_in_non_snp_event = findrecombinantSNPs(counts,ancnti,p,goodpos)

involved_in_non_snp_event=false(size(goodpos));
goodpos = find(goodpos);
Nsample = size(counts,3);

ancnti_m=repmat(ancnti,1,Nsample);
%
[cmajorAF, cmajorNT, cminorNT, cminorAF] = div_major_allele_freq(counts);
cminorAF(isnan(cminorAF))=0;
%
mutantAF=zeros(size(cmajorNT));
mutantAF(cmajorNT~=ancnti_m)=cmajorAF(cmajorNT~=ancnti_m);
mutantAF(cminorNT~=ancnti_m)=mutantAF(cminorNT~=ancnti_m)+cminorAF(cminorNT~=ancnti_m); %this construction allows for positions with two different mutations

nonsnp=[];
distance_for_nonsnp=500;
for j=1:numel(goodpos)
    gp=p(goodpos(j));
    %find nearby snps
    region=find(p>gp-distance_for_nonsnp&p<gp+distance_for_nonsnp);
    if numel(region)>1
        r = mutantAF(region,:);
        corrmatrix = corr(r') ;
        [a,b]=find(corrmatrix>.75);
        nonsnp=[nonsnp region(a(a~=b))'];
    end
end
%

nonsnp=unique(nonsnp);
involved_in_non_snp_event(nonsnp)=true;

end