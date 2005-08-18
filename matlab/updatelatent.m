function [sBar, Sigma_s] = updatelatent(model, X)

% UPDATELATENT Update the source values in the variational algorithm.

% GCA

sBar = zeros(model.numData, model.latentDim);
Sigma_s = zeros(model.latentDim, model.latentDim, model.numData);

B = diag(model.beta);
ATBA = model.A'*B*model.A;

for n = 1:model.numData
  invSigma_s = diag(model.tau(n, :)) + ATBA;
  C = chol(invSigma_s);
  Cinv = eye(model.latentDim)/C;
  Sigma_s(:, :, n) = Cinv*Cinv'; 
  sBar(n, :) = (X(n, :).*model.beta)*model.A*Sigma_s(:, :, n);
end

