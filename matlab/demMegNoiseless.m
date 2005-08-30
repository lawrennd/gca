% DEMMEGNOISELESS Run the noiseless algorithm on HUT MEG data.

% GCA

rand('seed', 3.14e5)
randn('seed', 3.14e5);

load 'datasets/hutmeg.mat'
expectedLatentDim = 6;


% Extract a 20 second portion 24 seconds in (sample rate is 148.5 Hz)
X = MEG_art(:, 3500:6500)';
% Centre the data
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