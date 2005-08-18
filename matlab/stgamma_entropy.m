function h = stgamma_entropy(nu, sigma2);

% STGAMMA_ENTROPY Entropy of a Student-t parameterised gamma distribution.

% GCA

a=nu/2;
b = a.*sigma2;

h = gamma_entropy(a, b);
