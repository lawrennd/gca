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
global FANOISE
global DATADIM
global LATENTDIM
global NDATA

DATADIM = size(X, 2);
NDATA = size(X, 1);
LATENTDIM = size(SBAR, 2);

% Work out parameters from expectations.
TAU = 1./SIGMA2BAR_TAU;
for j = 1:LATENTDIM
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
	  - log(NUBAR_TAU(j)/2) ...
	  - log(SIGMA2BAR_TAU(:,j));
end




% Now do the expectations of the joint distribution

% Data term
expectedOutput = (X - SBAR*A');
datapart = 0;
if FANOISE
  for i = 1:DATADIM
    datapart = datapart + BETA(i)*expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:NDATA
      datapart = datapart + BETA(i)*A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
    end
  end

  datapart = datapart*0.5;
  expLnp_x = -datapart  - 0.5*NDATA*DATADIM*log(2*pi) + 0.5*NDATA*sum(log(BETA));
else
  for i = 1:DATADIM
    datapart = datapart + BETA*expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:NDATA
      datapart = datapart + BETA*A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
    end
  end

  datapart = datapart*0.5;
  expLnp_x = -datapart  - 0.5*NDATA*DATADIM*log(2*pi) + 0.5*NDATA*DATADIM*log(BETA);
end  
% Latent variables' prior
latentPart = 0;
for n = 1:NDATA
  latentPart = latentPart + SBAR(n, :)*diag(TAU(n, :))*SBAR(n, :)' + ...
      trace(diag(TAU(n, :))*SIGMA_S(:, :, n));
end
expLnp_s = -0.5*latentPart - 0.5*NDATA*LATENTDIM*log(2*pi) + 0.5*sum(sum(LNTAU));
  
% Noise variance prior
if FANOISE
  expLnp_BETA = 0;%-sum(log(BETA));
else
  expLnp_BETA = 0;%-log(BETA);
end

% TAU prior
expLnp_TAU = 0;	  
for j = 1:LATENTDIM
  expLnp_TAU = expLnp_TAU + stgamma_expectation(NU_TAU(j), ...
						SIGMA2_TAU(j), ...
						TAU(:, j), ...
						LNTAU(:, j));
end

negenergy = expLnp_x + expLnp_s + expLnp_TAU + expLnp_BETA; 

entropy_S = 0;
for n = 1:NDATA
  entropy_S = entropy_S+0.5*(log(det(SIGMA_S(:, :, n))) ...
		      + 0.5*LATENTDIM*(1+log(2*pi)));
end


entropy_TAU = 0;	  
for j = 1:LATENTDIM
  entropy_TAU = entropy_TAU + stgamma_entropy(NUBAR_TAU(j), SIGMA2BAR_TAU(:, j));
end


entropy = entropy_S + entropy_TAU;


lll = negenergy + entropy;






