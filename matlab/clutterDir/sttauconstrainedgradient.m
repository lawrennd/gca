function g = sttauconstrainedgradient(param, expTau, expLnTau, min_tau)

% BICA_TAUOBJECTIVE Lower bound for a Bayesian ICA as a function of tau parameters
% 
%	Description
%       function lll = bica_objective(tau, model, data)
%	See Also

%	Copyright (c) Neil D Lawrence (2001)
% Compute the parameters of the distributions

ndata = size(expTau, 1);
latentDim = size(expTau, 2);

g = zeros(size(param));
gamma = param(1:latentDim);
nu_tau = min_tau+exp(gamma);
lnnu_tau = log(nu_tau);
lnsigma2_tau = param(latentDim+1:2*latentDim);
sigma2_tau = exp(lnsigma2_tau);
lambda = param(2*latentDim+1:end);
scaleToo = 1;
sigma2_tau(find(sigma2_tau<eps)) = eps;

a = nu_tau/2;

for k = 1:latentDim
  g(k) = ndata/2*(log(a(k)) + 1 - psi(a(k)) + log(sigma2_tau(k))) + ...
       0.5*(sum(expLnTau(:, k)) - sigma2_tau(k)* ...
	    sum(expTau(:, k)));
  g(k+latentDim) = ndata*a(k)/(sigma2_tau(k)) - a(k)* ...
      sum(expTau(:, k));
end


g = -g.*[exp(gamma) sigma2_tau zeros(1, latentDim)];

for k = 1:latentDim
  g(k) = g(k) - exp(gamma(k))*lambda(k)*(1/(nu_tau(k)-2) - nu_tau(k)/(nu_tau(k)-2)^2)*sigma2_tau(k);
  g(k+latentDim) = g(k+latentDim) - lambda(k)*(nu_tau(k)/(nu_tau(k)-2)*sigma2_tau(k));
%  g(k+2*latentDim) = -(2/(exp(gamma(k))+1)*sigma2_tau(k) + 1;
  g(k+2*latentDim) = -(nu_tau(k)/(nu_tau(k)-2)*sigma2_tau(k)-1);
end
