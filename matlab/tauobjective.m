function f = tauobjective(lnparam, expTau, expLnTau, sigma2_tau)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of log tau parameters
% 
%	Description
%       function lll = bica_objective(lnnu_tau, latentDim, sigma2_tau, sigma2bar_tau, nubar_tau)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions


latentDim = size(expTau, 2);
c = exp(lnparam(1:latentDim));
d = exp(lnparam(latentDim+1:end));

f = 0;
for k = 1:latentDim
  
  f = f - gamma_expectation(c(k), ...
			   d(k), ...
			   expTau(:, k), ...
			   expLnTau(:, k));
end


