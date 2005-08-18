function [sbar, Sigma_s] = updatelatent2(XU, expD, expD2, V, tau, beta) 

ndata = size(XU, 1);
latentDim = size(tau, 2);
sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);

for n = 1:ndata
  invSigma_s = diag(tau(n, :)) + beta*V*expD2*V';
  U = chol(invSigma_s);
  s = warning;
  warning('off');
  Uinv = eye(latentDim)/U;
  warning(s);
  Sigma_s(:, :, n) = Uinv*Uinv'; 
  sbar(n, :) = XU(n, :)*beta*expD*V'*Sigma_s(:, :, n);
end
