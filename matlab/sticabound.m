function lll = sticabound(V, lambda, XULinv, beta, nu, sigma2_tau, sbar, ...
			   Sigma_s,  nubar_tau, sigma2bar_tau)

dataDim = size(XULinv, 2);
ndata = size(XULinv, 1);
latentDim = size(sbar, 2);

% Work out parameters from expectations.
tau = 1./sigma2bar_tau;
for j = 1:latentDim
  lntau(:, j) = psi(nubar_tau(j)/2) ...
	  - log(nubar_tau(j)/2) ...
	  - log(sigma2bar_tau(:,j));
end




% Now do the expectations of the joint distribution

% Data term
expectedOutput = (XULinv - sbar*V);
datapart = 0;
for i = 1:dataDim
  datapart = datapart + beta*lambda(i)*lambda(i)*expectedOutput(:, i)'*expectedOutput(:, i);
  for n = 1:ndata
    datapart = datapart + beta*lambda(i)*lambda(i)*V(:, i)'*Sigma_s(:, :, n)*V(:, i);
  end
end

datapart = datapart*0.5;
expLnp_x = -datapart  - 0.5*ndata*dataDim*log(2*pi) + 0.5*ndata*dataDim*log(beta);

% Latent variables' prior
latentPart = 0;
for n = 1:ndata
  latentPart = latentPart + sbar(n, :)*diag(tau(n, :))*sbar(n, :)' + ...
      trace(diag(tau(n, :))*Sigma_s(:, :, n));
end
expLnp_s = -0.5*latentPart - 0.5*ndata*latentDim*log(2*pi) + 0.5*sum(sum(lntau));
  
% Noise variance prior
expLnp_beta = -log(beta);

% Tau prior
expLnp_tau = 0;	  
for j = 1:latentDim
  expLnp_tau = expLnp_tau + stgamma_expectation(nu(j), ...
						sigma2_tau(j), ...
						tau(:, j), ...
						lntau(:, j));
end

negenergy = expLnp_x + expLnp_s + expLnp_tau + expLnp_beta; 

entropy_S = 0;
for n = 1:ndata
  entropy_S = entropy_S+0.5*(log(det(Sigma_s(:, :, n))) ...
		      + 0.5*latentDim*(1+log(2*pi)));
end


entropy_tau = 0;	  
for j = 1:latentDim
  entropy_tau = entropy_tau + stgamma_entropy(nubar_tau(j), sigma2bar_tau(:, j));
end


entropy = entropy_S + entropy_tau;


lll = negenergy + entropy;






