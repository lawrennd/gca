function f = priorobjectivesq(sqrtparam, expTau, expLnTau, min_tau, k)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

latentDim = size(expTau, 2);

g = zeros(size(sqrtparam));
gamma = sqrtparam(1);
if gamma*gamma < eps;
  gamma = sqrt(eps)
end
nu_tau = min_tau+gamma*gamma;
sqrtsigma2_tau = sqrtparam(2);
sigma2_tau = sqrtsigma2_tau*sqrtsigma2_tau;
lambda = sqrtparam(3);


f = - stgamma_expectation(nu_tau, ...
			  sigma2_tau, ...
			  expTau(:, k), ...
			  expLnTau(:, k)) - lambda*(nu_tau/(nu_tau-2)*sigma2_tau-1);











