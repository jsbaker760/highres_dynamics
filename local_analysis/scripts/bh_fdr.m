function [sig_bool,crit_p]= bh_fdr(p_vals,alpha)

% resize p-vals to horizontal array if vertical
if size(p_vals,1)>1&size(p_vals,2)==1
    p_vals=p_vals';
end

M = numel(p_vals);

[p_sorted]= sort(p_vals,'ascend');

crit_p = max(p_sorted(p_sorted<(((1:M)/M)*alpha)));

if isempty(crit_p)
    sig_bool = false(1,M);
else
    sig_bool=p_vals<=crit_p;
end
