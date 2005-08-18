function model = variationalGCA(model, X, lltol, calcLLEvery, min_tau, ...
				display, ppcaInit, maxIter)

% VARIATIONALGCA Run the variational GCA algorithm.

% GCA

if nargin < 8
  maxIter = 1000;
end
if nargin < 7
  ppcaInit = 0;
end
if nargin < 6
  display = 0;
end
% Initialisation of A
if ppcaInit
  if model.latentDim < model.dataDim
    if display
      fprintf('Initialising A with PPCA.\n');
    end
    % If it's underdetermined use ppca
    [variance, U, V] = ppca(cov(X), model.latentDim);
    model.A = U*diag(sqrt(V));
  else
    error('Can only initialise with ppca for underdetermined systems.');
  end
else
  % Otherwise random little values
  if display
    fprintf('Initialising A with small random values.\n');
  end
  model.A = randn(model.dataDim, model.latentDim)*0.1;
  variance = mean(var(X));
end

% Initialisation of prior
if ppcaInit
  nuInit = 100;
else
  nuInit = 5;
end
if display
  fprintf('Initialising nu_tau as %4.2f.\n', nuInit);
  fprintf('Minimum value for nu_tau is %4.2f.\n', min_tau);
end
model.nu_tau = nuInit*ones(1, model.latentDim);
model.sigma2_tau = (model.nu_tau - 2)./model.nu_tau;

% Initialisation of beta
if model.FANoise
  model.beta = nrepmat(1./variance, 1, model.dataDim);
else
  model.beta = 1./variance;
end

% Initialisation of taum
model.nuBar_tau = nrepmat(model.nu_tau, 1, model.numData);
model.sigma2Bar_tau = nrepmat(model.sigma2_tau, 1, model.numData);

model.tau = 1./model.sigma2Bar_tau;
for j = 1:model.latentDim
   model.lntau(:, j) = digamma(model.nuBar_tau(j)/2) ...
      - log(model.nuBar_tau(j)/2) ...
      - log(model.sigma2Bar_tau(:,j));
end

% If displaying, store the values
if display > 1
  Astore = model.A(:)';
  betaStore = model.beta(:)';
  nu_tauStore = model.nu_tau(:)';
end

% Compute S from other initialisations
[model.sBar, model.Sigma_s] = updatelatent(model, X);

% Computer tau from the new S values
[model.sigma2Bar_tau, model.nuBar_tau] = stupdatetau(model, X);
model.tau = 1./model.sigma2Bar_tau;
for j = 1:model.latentDim
   model.lntau(:, j) = digamma(model.nuBar_tau(j)/2) ...
      - log(model.nuBar_tau(j)/2) ...
      - log(model.sigma2Bar_tau(:,j));
end

counter = 1;
lll(counter) = sticabound(model, X);
llldiff = 1;

iter = 0;
while  iter < maxIter
  order = randperm(5);
  iter = iter + 1;
  for i = order
    switch i
     case 1
      % Update the latent variables
      [model.sBar, model.Sigma_s] = updatelatent(model, X);
     case 2
      % Update the taus
      [model.sigma2Bar_tau, model.nuBar_tau] = stupdatetau(model, X);
      model.tau = 1./model.sigma2Bar_tau;
      for j = 1:model.latentDim
	model.lntau(:, j) = digamma(model.nuBar_tau(j)/2) ...
	    - log(model.nuBar_tau(j)/2) ...
	    - log(model.sigma2Bar_tau(:,j));
      end
     case 3
      % Update the mixing matrix
      model.A = updateA(model, X);
     case 4
      % Update the noise variance
      model.beta = updatebeta(model, X);
     case 5 
      % Update the latent variable's parameters
      [model.sigma2_tau, model.nu_tau] = stupdatetauprior(model, X, 'scg', min_tau);
    end
  end
  
  if display > 1
    Astore(iter, :) = model.A(:)';
    betaStore(iter, :) = model.beta(:)';
    nu_tauStore(iter, :) = model.nu_tau(:)';
  end
  
  if ~rem(iter, calcLLEvery)
    counter = counter + 1;
    lll(counter) = sticabound(model, X);
    llldiff = lll(counter) - lll(counter-1);
    if llldiff < lltol*calcLLEvery
      break
    end
    if display
      fprintf('Iteration %i.%i, log likelihood change: %d\n', iter, i,llldiff)
    end
    
  
    
    if display > 1
      figure(1)
      subplot(4, 1, 1)
      plot(Astore)
      subplot(4, 1, 2)
      plot(log10(betaStore))
      subplot(4, 1, 3)
      plot(log10(nu_tauStore))
      subplot(4, 1, 4)
      plot(0:calcLLEvery:iter, lll(1:counter));
      drawnow
    end 
  end
end
