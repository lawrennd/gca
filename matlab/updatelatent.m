function [sbar, Sigma_s] = updatelatent(XUL, V, tau, beta, lambda, choose) 

ndata = size(XUL, 1);
if nargin < 6
  choose = 1:ndata;
end
latentDim = size(tau, 2);
dataDim = size(V, 1);
sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);
B = diag(beta*lambda.*lambda);
VBVT = V*B*V';
for n = choose
  invSigma_s = diag(tau(n, :)) + VBVT;
  C = chol(invSigma_s);
  Cinv = eye(latentDim)/C;
  Sigma_s(:, :, n) = Cinv*Cinv'; 
  sbar(n, :) = XUL(n, :)*B*V'*Sigma_s(:, :, n);
end
