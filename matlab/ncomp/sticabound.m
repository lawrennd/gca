function lll = sticabound

global A
global X
global BETA
global NU_TAU
global SIGMA2_TAU
global SBAR
global SIGMA_S
global NUBAR_TAU
global SIGMA2BAR_TAU
global TAU
global LNTAU

dataDim = size(X, 2);
ndata = size(X, 1);
latentDim = size(SBAR, 2);

% Work out parameters from expectations.
TAU = 1./SIGMA2BAR_TAU;
for j = 1:latentDim
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
	  - log(NUBAR_TAU(j)/2) ...
	  - log(SIGMA2BAR_TAU(:,j));
end




% Now do the expectations of the joint distribution

% Data term
expectedOutput = (X - SBAR*A');
datapart = 0;
for i = 1:dataDim
  datapart = datapart + BETA*expectedOutput(:, i)'*expectedOutput(:, i);
  for n = 1:ndata
    datapart = datapart + BETA*A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
  end
end

datapart = datapart*0.5;
expLnp_x = -datapart  - 0.5*ndata*dataDim*log(2*pi) + 0.5*ndata*dataDim*log(BETA);

% Latent variables' prior
latentPart = 0;
for n = 1:ndata
  latentPart = latentPart + SBAR(n, :)*diag(TAU(n, :))*SBAR(n, :)' + ...
      trace(diag(TAU(n, :))*SIGMA_S(:, :, n));
end
expLnp_s = -0.5*latentPart - 0.5*ndata*latentDim*log(2*pi) + 0.5*sum(sum(LNTAU));
  
% Noise variance prior
expLnp_BETA = -log(BETA);

% TAU prior
expLnp_TAU = 0;	  
for j = 1:latentDim
  expLnp_TAU = expLnp_TAU + stgamma_expectation(NU_TAU(j), ...
						SIGMA2_TAU(j), ...
						TAU(:, j), ...
						LNTAU(:, j));
end

negenergy = expLnp_x + expLnp_s + expLnp_TAU + expLnp_BETA; 

entropy_S = 0;
for n = 1:ndata
  entropy_S = entropy_S+0.5*(log(det(SIGMA_S(:, :, n))) ...
		      + 0.5*latentDim*(1+log(2*pi)));
end


entropy_TAU = 0;	  
for j = 1:latentDim
  entropy_TAU = entropy_TAU + stgamma_entropy(NUBAR_TAU(j), SIGMA2BAR_TAU(:, j));
end


entropy = entropy_S + entropy_TAU;


lll = negenergy + entropy;






