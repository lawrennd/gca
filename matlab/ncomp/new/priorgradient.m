function g = priorgradient(lnparam, expTau, expLnTau, min_tau, k)

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
gamma = lnparam(1);
nu_tau = min_tau+exp(gamma);
lnsigma2_tau = lnparam(2);
sigma2_tau = exp(lnsigma2_tau);
lambda = lnparam(3);
a = nu_tau/2;
if any(sigma2_tau < eps)
  sigma2_tau(find(sigma2_tau < eps)) = eps;
end

g(1) = ndata/2*(log(a) + 1 - psi(a) + log(sigma2_tau)) + ...
       0.5*(sum(expLnTau(:, k)) - sigma2_tau* ...
	    sum(expTau(:, k)))  ...
       + lambda*(1/(nu_tau-2) ...
		    - nu_tau/(nu_tau-2)^2)*sigma2_tau;
g(1) = -g(1).*exp(gamma);
g(2) = ndata*a/(sigma2_tau) - a* ...
    sum(expTau(:, k))  + lambda*(nu_tau/(nu_tau-2));
g(2) = -g(2).*sigma2_tau;
g(3) = -nu_tau/(nu_tau-2)*sigma2_tau + 1;







