function lll = bicabound(V, Sigma_V, lambda, XULinv, c, d, a_beta, b_beta, a_alpha, b_alpha, sbar, ...
			   Sigma_s,  cbar, dbar, abar_beta, ...
			 bbar_beta, abar_alpha, bbar_alpha)

dataDim = size(XULinv, 2);
ndata = size(XULinv, 1);
latentDim = size(sbar, 2);

% Work out parameters from expectations.
for j = 1:latentDim
  tau(:, j) = cbar(j)./dbar(:, j);
  lntau(:, j) = psi(cbar(j)) - log(dbar(:, j));
end



beta = abar_beta./bbar_beta;
lnbeta = psi(abar_beta) - log(bbar_beta); 

alpha = abar_alpha./bbar_alpha;
lnalpha = psi(abar_alpha) - log(bbar_alpha); 

for i = 1:dataDim
  VVT(:, :, i) = V(:, i)*V(:, i)' + Sigma_V(:, :, i);
end

% Now do the expectations of the joint distribution



% Data term
expectedOutput = (XULinv - sbar*V);
datapart = 0;
for i = 1:dataDim
  datapart = datapart + beta*lambda(i)*lambda(i)*expectedOutput(:, i)'*expectedOutput(:, i);
  for n = 1:ndata
    datapart = datapart + beta*lambda(i)*lambda(i)*(V(:, i)'*Sigma_s(:, :, n)*V(:, i) ...
	   + trace(Sigma_s(:, :, n)*Sigma_V(:, :, i)) ...
	   + sbar(n, :)*Sigma_V(:, :, i)*sbar(n, :)');
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
  expLnp_tau = expLnp_tau + gamma_expectation(c(j), ...
					      d(j), ...
					      tau(:, j), ...
					      lntau(:, j));
end

% Mixing prior
expLnp_alpha = gamma_expectation(a_alpha, b_alpha, alpha, lnalpha);


expLnp_V = 0;
for j = 1:latentDim
  expLnp_V = expLnp_V - alpha(j)*sum(VVT(j, j, :), 3);
end
expLnp_V = dataDim/2*(sum(lnalpha) - latentDim*log(2*pi)) + 0.5*expLnp_V;


negenergy = expLnp_x + expLnp_s + expLnp_tau + expLnp_beta + expLnp_V; 

entropy_S = 0;
for n = 1:ndata
  entropy_S = entropy_S+0.5*(log(det(Sigma_s(:, :, n))) ...
		      + 0.5*latentDim*(1+log(2*pi)));
end

entropy_beta = gamma_entropy(abar_beta, bbar_beta);

entropy_alpha = gamma_entropy(abar_alpha, bbar_alpha);

entropy_tau = 0;	  
for j = 1:latentDim
  entropy_tau = entropy_tau + gamma_entropy(cbar(j), dbar(:, j));
end

entropy_V = 0;
for i = 1:dataDim
  entropy_V = entropy_V + 0.5*log(det(Sigma_V(:, :, i)));
end
entropy_V = entropy_V + dataDim*latentDim*(1+log(2*pi));


entropy = entropy_S + entropy_beta + entropy_tau + entropy_V + entropy_alpha;


lll = negenergy + entropy;






