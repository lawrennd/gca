function f = priorobjective(lnparam, expTau, expLnTau, min_tau, k)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

latentDim = size(expTau, 2);

g = zeros(size(lnparam));
gamma = lnparam(1);
nu_tau = min_tau+exp(gamma);
lnsigma2_tau = lnparam(2);
sigma2_tau = exp(lnsigma2_tau);
lambda = lnparam(3);
if any(sigma2_tau < eps)
  sigma2_tau(find(sigma2_tau < eps)) = eps;
end


f = 0;
f = f - stgamma_expectation(nu_tau, ...
			    sigma2_tau, ...
			    expTau(:, k), ...
			    expLnTau(:, k)) - lambda*(nu_tau/(nu_tau-2)*sigma2_tau-1);


