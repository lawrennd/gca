function d = updatetaupriormoments(tau, c)

a = 0;
b = 0;
ndata = size(tau, 1);
d = (a+ndata*c)./(b+sum(tau));

