function Sigma_s = updateSigma_s(V, tau, beta, lambda, choose) 

ndata = size(tau, 1);
if nargin < 5
  choose = 1:ndata;
end
latentDim = size(tau, 2);
dataDim = size(V, 1);
Sigma_s = zeros(latentDim, latentDim, ndata);
B = diag(beta*lambda.*lambda);
VBVT = V*B*V';
for n = choose
  invSigma_s = diag(tau(n, :)) + VBVT;
  C = chol(invSigma_s);
  Cinv = eye(latentDim)/C;
  Sigma_s(:, :, n) = Cinv*Cinv'; 
end
