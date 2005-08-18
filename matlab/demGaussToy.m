% DEMGAUSSTOY run the variational algorithm on the Gaussian toy problem.

% GCA

randn('seed', 1e5);
rand('seed', 1e5);

ndata = 1000;
varVector = 10:-1:1;
dataDim = length(varVector);

% Choose a projection
tempMat = randn(dataDim, 10);
[U, V] = eig(tempMat*tempMat');
project = U;

% Sample a Gaussian data set
gaussX = randn(ndata, dataDim);
for i = 1:dataDim;
  gaussX(:, i) = gaussX(:, i)*sqrt(varVector(i));
end
gaussX = gaussX*project';

global model.numData
global model.dataDim
global model.latentDim

global X

global model.beta
global model.nu_tau
global model.sigma2_tau

global model.nuBar_tau
global model.sigma2Bar_tau

global model.tau
global  model.lntau

global model.A
global model.sBar
global model.Sigma_s

global Fmodel.ANOISE % Set non-zero to use factor analysis noise model

display = 2;
X = gaussX;

Fmodel.ANOISE = 0;
model.latentDim = 5;

model.numData = size(X, 1);
model.dataDim = size(X, 2);

% X = zero meaned data;
meanX = mean(X, 1);
for i = 1:model.numData
  X(i, :)  = X(i, :) - meanX;
end


% Tolerances
lltol = 1e-5;     % The tolerance on the log-likelihood for convergence
calcLLEvery = 25; % How often to evaluate bound on log-likelihood
min_tau = 2.001; % The minimum value allowed for model.nu_tau

% Initialisations
model.A = randn( model.dataDim, model.latentDim)*0.01;
model.nu_tau = 5*ones(1, model.latentDim);
model.sigma2_tau = (model.nu_tau - 2)./model.nu_tau;

if Fmodel.ANOISE
  model.beta = 1./var(X);
else
  model.beta = model.dataDim/sum(var(X));
end

model.nuBar_tau = nrepmat(model.nu_tau, 1, model.numData);
model.sigma2Bar_tau = nrepmat(model.sigma2_tau, 1, model.numData);

model.tau = 1./model.sigma2Bar_tau;
for j = 1:model.latentDim
   model.lntau(:, j) = psi(model.nuBar_tau(j)/2) ...
      - log(model.nuBar_tau(j)/2) ...
      - log(model.sigma2Bar_tau(:,j));
end

if display > 1
  model.Astore = model.A(:)';
  betaStore = model.beta(:)';
  nu_tauStore = model.nu_tau(:)';
end

model.sBar = zeros(model.numData, model.latentDim);
model.Sigma_s = zeros(model.latentDim, model.latentDim, model.numData);
updatelatent(model, X)

stupdatetau(model, X)
model.tau = 1./model.sigma2Bar_tau;
for j = 1:model.latentDim
   model.lntau(:, j) = psi(model.nuBar_tau(j)/2) ...
      - log(model.nuBar_tau(j)/2) ...
      - log(model.sigma2Bar_tau(:,j));
end

counter = 1;
lll(counter) = sticabound(model, X)/model.numData;
llldiff = 1;

iter = 0;
while  iter < 10000;
  order = randperm(5);
  iter = iter + 1;
  for i = order
    switch i
     case 1
      updatelatent(model, X)
     case 2
      stupdatetau(model, X)
      model.tau = 1./model.sigma2Bar_tau;
      for j = 1:model.latentDim
	 model.lntau(:, j) = psi(model.nuBar_tau(j)/2) ...
	    - log(model.nuBar_tau(j)/2) ...
	    - log(model.sigma2Bar_tau(:,j));
      end
     case 3
      updatemodel.A(model, X)
     case 4
      updatebeta(model, X)
     case 5 
      stupdatetau(model, X)prior('scg', min_tau)
    end
  end

  
  if ~rem(iter, calcLLEvery)
    counter = counter + 1;
    lll(counter) = sticabound(model, X)/model.numData;
    llldiff = lll(counter) - lll(counter-1);
    if llldiff < lltol*calcLLEvery
      break
    end
    if display
      fprintf('Iteration %i.%i, log likelihood change: %d\n', iter, i,llldiff)
    end
  end

  
  if display > 1
    model.Astore(iter, :) = model.A(:)';
    betaStore(iter, :) = model.beta(:)';
    nu_tauStore(iter, :) = model.nu_tau(:)';
  end
  
  if display > 1
    figure(1)
    subplot(4, 1, 1)
    plot(model.Astore)
    subplot(4, 1, 2)
    plot(log10(betaStore))
    subplot(4, 1, 3)
    plot(log10(nu_tauStore))
    subplot(4, 1, 4)
    plot(1:counter, lll(1:counter));
    drawnow
  end  
end



model.Abar = model.A./nrepmat(sqrt(sum(model.A.*model.A)), 1, 10);
P = model.Abar'*project;
projections = sqrt(sum(P.^2));


% PPCmodel.A experiment
[var, Upca, lambda] = ppca(cov(X), 5);
Ppca = Upca'*project;
projectionsPPCmodel.A = sqrt(sum(Ppca.^2));

[U I] = eig(cov(X));
PCmodel.A = U'*project;
projectionsPCmodel.A = sqrt(sum(PCmodel.A.^2));

