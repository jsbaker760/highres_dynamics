function D = calc_alpha_div(D,m)
m = string(m);
em = sum(D)==0;
if em
    D=NaN;
elseif m=="shannon"
    D = -1*nansum(D.*log(D));
elseif m=="simpson"
    D = 1-(sum(D.^2));
end
