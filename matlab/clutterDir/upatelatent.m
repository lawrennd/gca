function [sbar, Sigma_s] = upatelatent(XUL, V, tau, beta, lambda) 

ndata = size(XUL, 1);
latentDim = size(tau, 2);
dataDim = size(V, 1);
sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);
B = diag(beta*lambda.*lambda);
for n = 1:ndata
  
  invSigma_s = diag(tau(n, :)) + V*B*V';
  C = chol(invSigma_s);
  Cinv = eye(latentDim)/C;
  Sigma_s(:, :, n) = Cinv*Cinv'; 
  sbar(n, :) = XUL(n, :)*B*V'*Sigma_s(:, :, n);
end
