function model = noiselessGCA(model, X, lltol, min_tau, display, pcaInit);

% NOISELESSGCA Run the noiseless GCA algorithm.

% GCA

if nargin < 6
  pcaInit = 0;
end
if nargin < 5
  display = 0;
end

% Initialisation of A
if pcaInit
  if display
    fprintf('Initialising inv(A) with PCA.\n');
  end
  % If it's underdetermined use ppca
  [U, V] = eig(cov(X));
  model.nu_tau = 5*ones(1, model.latentDim);
  model.sigma2_tau = (model.nu_tau-2)./model.nu_tau;
  W = U*sqrt(diag(1./diag(V)));
  W = W';
else
  % Otherwise random little values
  if display
    fprintf('Initialising inv(A) with small random values.\n');
  end
  W = randn(model.latentDim, model.dataDim)*0.1;
end

% Initialisation of nu_tau
% Initialisation of prior
if pcaInit
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


options = foptions;
options(1) = display;
options(3) = lltol;
options(14) = 10000;

% pack the parameters into a vector
params = W(:)';
params = [params sqrt(model.nu_tau - min_tau)];

% optimise using non-linear optimiser
params = quasinew('noiselessobjective', params, options, 'noiselessgradient', X, min_tau);

% unpack the weights
W = reshape(params(1:model.latentDim*model.latentDim), model.latentDim, size(X, 2));
model.sBar = X*W';

% Compute the mixing matrix
model.A = inv(W);

% Unpack the nu_taus
model.nu_tau = min_tau + params(model.latentDim*model.latentDim+1:end).^2;
model.sigma2_tau = (model.nu_tau-2)./model.nu_tau;