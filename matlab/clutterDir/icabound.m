function lll = icabound(V, lambda, XULinv, c, d, a_beta, b_beta, sbar, ...
			   Sigma_s,  cbar, dbar, abar_beta, bbar_beta)

dataDim = size(XULinv, 2);
ndata = size(XULinv, 1);
latentDim = size(sbar, 2);

% Work out parameters from expectations.
tau = 1./dbar;
for j = 1:latentDim
  lntau(:, j) = psi(cbar(j)/2) ...
	  - log(cbar(j)/2) ...
	  - log(dbar(:,j));
end



beta = abar_beta./bbar_beta;
lnbeta = psi(abar_beta) - log(bbar_beta); 

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
expLnp_x = -datapart  - 0.5*ndata*dataDim*log(2*pi) + 0.5*ndata*dataDim*sum(lnbeta);

% Latent variables' prior
latentPart = 0;
for n = 1:ndata
  latentPart = latentPart + sbar(n, :)*diag(tau(n, :))*sbar(n, :)' + ...
      trace(diag(tau(n, :))*Sigma_s(:, :, n));
end
expLnp_s = -0.5*latentPart - 0.5*ndata*latentDim*log(2*pi) + 0.5*sum(sum(lntau));
  
% Noise variance prior
expLnp_beta = gamma_expectation(a_beta, b_beta, beta, lnbeta);

% Tau prior
expLnp_tau = 0;	  
for j = 1:latentDim
  expLnp_tau = expLnp_tau + stgamma_expectation(c(j), ...
						d(j), ...
						tau(:, j), ...
						lntau(:, j));
end

negenergy = expLnp_x + expLnp_s + expLnp_tau + expLnp_beta; 

entropy_S = 0;
for n = 1:ndata
  entropy_S = entropy_S+0.5*(log(det(Sigma_s(:, :, n))) ...
		      + 0.5*latentDim*(1+log(2*pi)));
end

entropy_beta = gamma_entropy(abar_beta, bbar_beta);

entropy_tau = 0;	  
for j = 1:latentDim
  entropy_tau = entropy_tau + stgamma_entropy(cbar(j), dbar(:, ...
						  j));
end


entropy = entropy_S + entropy_beta + entropy_tau;


lll = negenergy + entropy;






