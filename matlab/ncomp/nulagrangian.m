function f = nulagrangian(nu, min_tau, j)

global NDATA

global SIGMA2_TAU
global LAGRANGE

global TAU
global LNTAU

%nu = min_tau + exp(gamma);
f = 0.5*NDATA*nu*log(.5*SIGMA2_TAU(j)*nu) ...
    - NDATA*gammaln(nu/2) ...
    + (nu/2 - 1)*sum(LNTAU(:, j)) ...
    - (nu/2)*SIGMA2_TAU(j)*sum(TAU(:, j)) ...
    + LAGRANGE(j)*(nu/(nu-2)*SIGMA2_TAU(j) - 1);
f=-f;