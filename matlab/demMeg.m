% DEMMEG Run the variational algorithm on HUT MEG data.

% GCA

rand('seed', 3.14e5)
randn('seed', 3.14e5);

load 'datasets/hutmeg.mat'
expectedLatentDim = 6;


% Centre the data
X = MEG_art(:, 3500:6500)';
meanX = mean(X, 1);
model.numData = size(X, 1);
model.dataDim = size(X, 2);
for i = 1:model.numData
  X(i, :)  = X(i, :) - meanX;
end

% Model options
model.FANoise = 0;
model.latentDim = 7;
min_tau = 2.5; % The minimum value allowed for model.nu_tau

% Algorithmic options
lltol = 1e-4;     % The tolerance on the log-likelihood for convergence
calcLLEvery = 100; % How often to evaluate bound on log-likelihood
display = 2;      % Display the results or not. 
maxIter = 2e4; % Maximum number of iterations
ppcaInit = 0; % initialise with ppca?

model = variationalGCA(model, X, lltol, calcLLEvery, min_tau, display, ...
		       ppcaInit, maxIter);