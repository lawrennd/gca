function lll = sticabound(model, X)

% STICABOUND Compute the lower bound on the log likelihood.

% GCA

model.dataDim = size(X, 2);
model.numData = size(X, 1);
model.latentDim = size(model.sBar, 2);

% Work out parameters from expectations.
model.tau = 1./model.sigma2Bar_tau;
for j = 1:model.latentDim
   model.lntau(:, j) = digamma(model.nuBar_tau(j)/2) ...
	  - log(model.nuBar_tau(j)/2) ...
	  - log(model.sigma2Bar_tau(:,j));
end


% Now do the expectations of the joint distribution

% Data term
expectedOutput = (X - model.sBar*model.A');
datapart = 0;
if model.FANoise
  for i = 1:model.dataDim
    datapart = datapart + model.beta(i)*expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:model.numData
      datapart = datapart + model.beta(i)*model.A(i, :)*model.Sigma_s(:, :, n)*model.A(i, :)';
    end
  end

  datapart = datapart*0.5;
  expLnp_x = -datapart  - 0.5*model.numData*model.dataDim*log(2*pi) + 0.5*model.numData*sum(log(model.beta));
else
  for i = 1:model.dataDim
    datapart = datapart + model.beta*expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:model.numData
      datapart = datapart + model.beta*model.A(i, :)*model.Sigma_s(:, :, n)*model.A(i, :)';
    end
  end

  datapart = datapart*0.5;
  expLnp_x = -datapart  - 0.5*model.numData*model.dataDim*log(2*pi) + 0.5*model.numData*model.dataDim*log(model.beta);
end  
% Latent variables' prior
latentPart = 0;
for n = 1:model.numData
  latentPart = latentPart + model.sBar(n, :)*diag(model.tau(n, :))*model.sBar(n, :)' + ...
      trace(diag(model.tau(n, :))*model.Sigma_s(:, :, n));
end
expLnp_s = -0.5*latentPart - 0.5*model.numData*model.latentDim*log(2*pi) + 0.5*sum(sum( model.lntau));
  
% Noise variance prior
if model.FANoise
  expLnp_beta = 0;%-sum(log(model.beta));
else
  expLnp_beta = 0;%-log(model.beta);
end

% model.tau prior
expLnp_tau = 0;	  
for j = 1:model.latentDim
  expLnp_tau = expLnp_tau + stgamma_expectation(model.nu_tau(j), ...
						model.sigma2_tau(j), ...
						model.tau(:, j), ...
						 model.lntau(:, j));
end

negenergy = expLnp_x + expLnp_s + expLnp_tau + expLnp_beta; 

entropy_S = 0;
for n = 1:model.numData
  entropy_S = entropy_S+0.5*(log(det(model.Sigma_s(:, :, n))) ...
		      + 0.5*model.latentDim*(1+log(2*pi)));
end


entropy_tau = 0;	  
for j = 1:model.latentDim
  entropy_tau = entropy_tau + stgamma_entropy(model.nuBar_tau(j), model.sigma2Bar_tau(:, j));
end


entropy = entropy_S + entropy_tau;


lll = negenergy + entropy;






