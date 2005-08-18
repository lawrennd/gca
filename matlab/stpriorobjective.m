function f = stpriorobjective(gamma, expTau, expLnTau, min_tau, k)

% STPRIOROBJECTIVE The likelihood lower bound as a function of gamma.

% GCA

latentDim = size(expTau, 2);

if gamma*gamma < eps;
  gamma = sqrt(eps);
end
nu_tau = min_tau+gamma*gamma;
sigma2_tau = (nu_tau-2)/nu_tau;


f = - stgamma_expectation(nu_tau, ...
			  sigma2_tau, ...
			  expTau(:, k), ...
			  expLnTau(:, k));











