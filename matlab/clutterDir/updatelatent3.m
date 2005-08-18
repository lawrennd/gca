function [sbar, Sigma_s] = updatelatent3(UX, expD, expD2, V, tau, beta) 

ndata = size(UX, 1);
latentDim = size(tau, 2);
sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);

for n = 1:ndata
  invSigma_s = diag(tau(n, :)) + V*expD2*diag(beta)*V';
  U = chol(invSigma_s);
  s = warning;
  warning('off');
  Uinv = eye(latentDim)/U;
  warning(s);
  Sigma_s(:, :, n) = Uinv*Uinv'; 
  sbar(n, :) = UX(n, :)*diag(beta)*expD*V'*Sigma_s(:, :, n);
end
