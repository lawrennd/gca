function [sbar, Sigma_s] = bayesupdatelatent(XULinv, V, Sigma_V, tau, beta, lambda) 

ndata = size(XULinv, 1);
latentDim = size(tau, 2);
dataDim = size(XULinv, 2);
sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);
B = diag(beta*lambda.*lambda);
for n = 1:ndata
  invSigma_s = diag(tau(n, :));
  for i = 1:dataDim
    invSigma_s = invSigma_s + B(i, i)*(V(:, i)*V(:, i)' + ...
	Sigma_V(:, :, i));
  end
  C = chol(invSigma_s);
  %s = warning;
  %warning('off');
  Cinv = eye(latentDim)/C;
  %warning(s);
  Sigma_s(:, :, n) = Cinv*Cinv'; 
  sbar(n, :) = XULinv(n, :)*B*V'*Sigma_s(:, :, n);
end
