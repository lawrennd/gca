% DEMFOETALNOISELESS Run the noiseless algorithm on the foetal ecg data.

% GCA 

rand('seed', 3.14e5)
randn('seed', 3.14e5);

load 'datasets/foetal_ecg.dat'

% Centre the data
X = foetal_ecg(:, 2:9);
meanX = mean(X, 1);
model.numData = size(X, 1);
model.dataDim = size(X, 2);
for i = 1:model.numData
  X(i, :)  = X(i, :) - meanX;
end

% Model options
model.dataDim = size(X, 2);
model.latentDim = model.dataDim;
model.numData = size(X, 1);
min_tau = 2.5;

% Algorithmic options
lltol = 1e-4;     % The tolerance on the log-likelihood for convergence
display = 1;      % Display the results or not. 
pcaInit = 0;      % Initialise with pca?

model = noiselessGCA(model, X, lltol, min_tau, display, pcaInit);