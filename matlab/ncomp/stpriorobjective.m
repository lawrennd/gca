function f = stpriorobjective(gamma, expTau, expLnTau, min_tau, k)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

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











