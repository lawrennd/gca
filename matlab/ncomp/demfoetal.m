%clear all;
HOME = getenv('HOME');
addpath([HOME '/mlprojects/gca/matlab'])
rand('seed', 3.14e5)
randn('seed', 3.14e5);
display = 2;
load 'c:\datasets\ICA\ECG Data\foetal_ecg.dat'

global NDATA
global DATADIM
global LATENTDIM

global X

global BETA
global NU_TAU
global NU_GAUSS
global SIGMA2_TAU

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global A
global SBAR
global SIGMA_S

global FANOISE % Set non-zero to use factor analysis noise model

X = foetal_ecg(:, 2:9);

FANOISE = 0;
LATENTDIM = 16;

NDATA = size(X, 1);
DATADIM = size(X, 2);

% X = zero meaned data;
meanX = mean(X, 1);
for i = 1:NDATA
  X(i, :)  = X(i, :) - meanX;
end


% Tolerances
lltol = 1e-5;     % The tolerance on the log-likelihood for convergence
calcLLEvery = 25; % How often to evaluate bound on log-likelihood
min_tau = 2.5; % The minimum value allowed for NU_TAU

% Initialisations
if LATENTDIM < DATADIM
    [variance, U, V] = ppca(cov(X), LATENTDIM);
    A = U*diag(sqrt(V));
else
    A = randn(DATADIM, LATENTDIM)*0.1;
    variance = mean(var(X));
end

NU_TAU = 5*ones(1, LATENTDIM);
SIGMA2_TAU = (NU_TAU - 2)./NU_TAU;

if FANOISE
  BETA = 1./variance;
else
  BETA = DATADIM/sum(var(X));
end

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
lll(counter) = sticabound/NDATA;
llldiff = 1;

iter = 0;
while  iter < 10000
  order = randperm(5);
  iter = iter + 1;
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

  if display > 1
    Astore(iter, :) = A(:)';
    betaStore(iter, :) = BETA(:)';
    nu_tauStore(iter, :) = NU_TAU(:)';
  end
  
  if ~rem(iter, calcLLEvery)
    counter = counter + 1;
    lll(counter) = sticabound/NDATA;
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
