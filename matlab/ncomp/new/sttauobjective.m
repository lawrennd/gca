function f = sttauobjective(lnparam, expTau, expLnTau, sigma2_tau, lambda, min_tau)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions


latentDim = size(expTau, 2);
if length(lnparam) >2*latentDim
  lambda = lnparam(latentDim*2+1:end);
end
if length(lnparam) > latentDim
  nu_tau = min_tau+exp(lnparam(1:latentDim));
  sigma2_tau = exp(lnparam(latentDim+1:latentDim*2));
else
  nu_tau = min_tau+exp(lnparam);
end

f = 0;
for j = 1:latentDim
  f = f - stgamma_expectation(nu_tau(j), ...
			   sigma2_tau(j), ...
			   expTau(:, j), ...
			   expLnTau(:, j)) - lambda(j)*(nu_tau(j)/(nu_tau(j)-2)*sigma2_tau(j)-1);
end


