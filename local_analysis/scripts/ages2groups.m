function groups=ages2groups(ages)

groups = zeros(size(ages));

groups(ages<=12) = 1;
groups(ages>12&ages<18) = 2;
groups(ages>29) = 3;

