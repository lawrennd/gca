function g = sttaugradient(lnparam, expTau, expLnTau, sigma2_tau, lambda, ...
			   min_tau)

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
if length(lnparam) > latentDim
  gamma = lnparam(1:latentDim);
  nu_tau = min_tau+exp(gamma);
  lnsigma2_tau = lnparam(latentDim+1:end);
  sigma2_tau = exp(lnsigma2_tau);
  
  scaleToo = 1;
else
  scaleToo = 0;
  gamma = lnparam;
  nu_tau = min_tau+exp(gamma);
  sigma2_tau = sigma2_tau;
end
a = nu_tau/2;
if any(sigma2_tau < eps)
  sigma2_tau(find(sigma2_tau < eps)) = eps;
end

for k = 1:latentDim
  g(k) = ndata/2*(log(a(k)) + 1 - psi(a(k)) + log(sigma2_tau(k))) + ...
       0.5*(sum(expLnTau(:, k)) - sigma2_tau(k)* ...
	    sum(expTau(:, k)))  ...
	 + lambda(k)*(1/(nu_tau(k)-2) ...
		      - nu_tau(k)/(nu_tau(k)-2)^2)*sigma2_tau(k);
  if scaleToo
    g(k+latentDim) = ndata*a(k)/(sigma2_tau(k)) - a(k)* ...
	sum(expTau(:, k))  + lambda(k)*(nu_tau(k)/(nu_tau(k)-2));
  end
end


if ~scaleToo
  g = -g.*exp(gamma);
else
  g = -g.*[exp(gamma) sigma2_tau];
end









