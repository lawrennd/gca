function g = taugradient(lnparam, expTau, expLnTau, sigma2_tau)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of tau parameters
% 
%	Description
%       function lll = bica_objective(tau, model, data)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

ndata = size(expTau, 1);
latentDim = size(expTau, 2);

g = zeros(size(lnparam));
c = exp(lnparam(1:latentDim));
d = exp(lnparam(latentDim+1:end));
  
for k = 1:latentDim
  g(k) = ndata*(log(d(k)) - psi(c(k))) + sum(expLnTau(:, k));
  g(k+latentDim) = ndata*c(k)/d(k) - sum(expTau(:, k));
end


g = -g.*[c d];




