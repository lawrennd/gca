%clear all;
HOME = getenv('HOME');
addpath([HOME '/mlprojects/gca/matlab'])
rand('seed', 1e2)
randn('seed', 1e2);
display = 1;
global X
load 'c:\datasets\faces\neilfaceData.mat'

global NDATA
global DATADIM
global LATENTDIM


global BETA
global NU_TAU
global SIGMA2_TAU

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global A
global SBAR
global SIGMA_S

global FANOISE % Set non-zero to use factor analysis noise model
X = double(X);

FANOISE = 0;
LATENTDIM = 100;

NDATA = size(X, 1);
DATADIM = size(X, 2);

% X = zero meaned data;
meanX = mean(X, 1);
for n = 1:NDATA
  X(n, :)  = X(n, :) - meanX;
end


% Tolerances
Atol = 1e-5;
lltol = 1e-5;     % The tolerance on the log-likelihood for convergence
calcLLEvery = 50; % How often to evaluate bound on log-likelihood
min_tau = 2.001; % The minimum value allowed for NU_TAU
[variance, U, lambda] = ppca(cov(X), LATENTDIM);
% Initialisations

A = U*diag(lambda);
NU_TAU = 5*ones(1, LATENTDIM);
SIGMA2_TAU = (NU_TAU - 2)./NU_TAU;

BETA = 1/variance;

NUBAR_TAU = nrepmat(NU_TAU, 1, NDATA);
SIGMA2BAR_TAU = nrepmat(SIGMA2_TAU, 1, NDATA);

TAU = 1./SIGMA2BAR_TAU;
for j = 1:LATENTDIM
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
      - log(NUBAR_TAU(j)/2) ...
      - log(SIGMA2BAR_TAU(:,j));
end

if display > 1
  Astore = A(:)';
  betaStore = BETA(:)';
  nu_tauStore = NU_TAU(:)';
end

SBAR = zeros(NDATA, LATENTDIM);
SIGMA_S = zeros(LATENTDIM, LATENTDIM, NDATA);
updatelatent

stupdatetau
TAU = 1./SIGMA2BAR_TAU;
for j = 1:LATENTDIM
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
      - log(NUBAR_TAU(j)/2) ...
      - log(SIGMA2BAR_TAU(:,j));
end

counter = 1;
%lll(counter) = sticabound/NDATA;
llldiff = 1;
oldA = A;
iter = 0;
lastorder = randperm(5);
while  iter < 10000
  iter = iter + 1;
  order = randperm(5);
  while(order(1) == lastorder(end))
    order = randperm(5);
  end
  lastorder = order;
  for i = order
    switch i
     case 1
      updatelatent
     case 2
      stupdatetau
      TAU = 1./SIGMA2BAR_TAU;
      for j = 1:LATENTDIM
	LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
	    - log(NUBAR_TAU(j)/2) ...
	    - log(SIGMA2BAR_TAU(:,j));
      end
     case 3
      updateA
     case 4
      updatebeta
     case 5 
      stupdatetauprior('scg', min_tau)
    end
  end
  Achange = max(max(abs(A-oldA)));
  oldA = A;
%  if ~rem(iter, calcLLEvery)
    counter = counter + 1;
%    lll(counter) = sticabound/NDATA;
%    llldiff = lll(counter) - lll(counter-1);
%    if Achange < Atol
%      break
%    end
%    if display
      fprintf('Iteration %i.%i, matrix change: %d\n', iter, i, Achange)
%    end
%  end

  
  if display > 1
    Astore(iter, :) = A(:)';
    betaStore(iter, :) = BETA(:)';
    nu_tauStore(iter, :) = NU_TAU(:)';
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
    plot(1:counter, lll(1:counter));
    drawnow
  end  
  drawnow
  if ~rem(iter, 10)
    save facesA A BETA NU_TAU
  end
end
