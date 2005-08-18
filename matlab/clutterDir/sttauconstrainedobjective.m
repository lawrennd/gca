function f = sttauconstrainedobjective(param, expTau, expLnTau, min_tau)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions


latentDim = size(expTau, 2);
nu_tau = min_tau+exp(param(1:latentDim));
sigma2_tau = exp(param(latentDim+1:2*latentDim));

lambda = param(2*latentDim+1:end);
f = 0;
for j = 1:latentDim
  f = f - stgamma_expectation(nu_tau(j), ...
			   sigma2_tau(j), ...
			   expTau(:, j), ...
			   expLnTau(:, j)) ...
  -lambda(j)*(nu_tau(j)/(nu_tau(j)-2)*sigma2_tau(j)-1);
end


