function L = icalikelihood(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
			   Sigma_s,  tau, lntau, beta, lnbeta)

dataDim = size(XUL, 2);
ndata = size(XUL, 1);
latentDim = size(sbar, 2);

% FIrst we do the expectations of the likelihood.

% Data term

expectedOutput = (XUL - sbar*V);
bbar_beta = zeros(1, dataDim);
for i = 1:dataDim
  bbar_beta(i) = bbar_beta(i) + expectedOutput(:, i)'*expectedOutput(:, i);
  for n = 1:ndata
    bbar_beta(i) = bbar_beta(i) + V(:, i)'*Sigma_s(:, :, n)*V(:, i);
  end
end
bbar_beta = 0.5*bbar_beta;
datapart = sum(bbar_beta);
datapart = datapart  - 0.5*ndata*dataDim*log(2*pi) + 0.5*ndata*sum(lnbeta);
bbar_beta = bbar_beta + b_beta;

% Latent variables' prior
latentPart = sbar
for n = 1:ndata
  latentPart = latentPart + sbar(n, :)*diag(tau)*sbar(n, :)' + ...
      trace(diag(tau)*Sigma_s(:, :, n));
end
latentPart = 0.5*latentPart - 0.5*ndata*latentDim*log(2*pi) + 0.5*sum(sum(lntau));
  
% Noise variance prior
expBetaPart = gamma_expectation(a_beta, b_beta, beta, lnbeta);

% Tau prior
expTauPart = gamma_expectation(nu, sigma2_tau, tau, lntau);

negenergy = expLnp_t + expLnp_x + expLnp_tau + expLnp_W + expLnp_alpha + ...
    expLnp_mu + expLnp_beta; 

entropy_S = 0;
for n = 1:ndata
  entropy_S = entropy_S+0.5*(log(det(Sigma_s(:, :, n))) ...
		      + 0.5*latentDim*(1+log(2*pi)));
end

abar_beta = a_beta + 0.5*ndata;

entropy_beta = gamma_entropy(abar_beta, bbar_beta);


entropy_tau = stgamma_entropy(nubar_tau, sigma2bar_tau);

entropy = entropy_S + entropy_W + entropy_beta + entropy_alpha + entropy_tau;
lll = negenergy + entropy;





% Next we do the entropies of the distributions

