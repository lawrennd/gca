function g = stpriorgradient(gamma, expTau, expLnTau, min_tau, k)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of tau parameters
% 
%	Description
%       function lll = bica_objective(tau, model, data)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

ndata = size(expTau, 1);
latentDim = size(expTau, 2);

if gamma*gamma < eps
  gamma = sqrt(eps);
end

nu_tau = min_tau+gamma*gamma;
a = nu_tau/2;

g(1) = ndata/2*(log(nu_tau-2) -log(2) -psi(a) + (nu_tau)/(nu_tau-2)) ...
       + 0.5*(sum(expLnTau(:, k)) - sum(expTau(:, k)));

g(1) = -g(1).*2*gamma;





