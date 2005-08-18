function exp = stgamma_expectation(nu, sigma2, expX, expLnX);

% STGAMMA_EXPECTATION Expectation of the Student-t parameterised gamma.

% GCA

a = nu/2;
b = a.*sigma2;

exp = gamma_expectation(a, b, expX, expLnX);
