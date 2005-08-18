function [V, Sigma_V] = bayesupdateV(sbar, Sigma_s, XULinv, beta, alpha, lambda)


dataDim = size(XULinv, 2);
latentDim = size(sbar, 2);
ndata = size(XULinv, 1);

sum_ssT = sbar'*sbar + sum(Sigma_s, 3);


for i = 1:dataDim
  invSigma_V = diag(alpha);
  invSigma_V = invSigma_V + beta*lambda(i)*lambda(i)*(sum_ssT);
  
  U = chol(invSigma_V);
  Uinv =  eye(latentDim)/U;
  Sigma_V(:, :, i) = Uinv*Uinv';
  V(:, i) = zeros(latentDim, 1);
  for n = 1:ndata
      V(:, i) = V(:, i) ...
	  + sbar(n, :)'*XULinv(n, i);
  end
  V(:, i) = Sigma_V(:, :, i)*beta*lambda(i)*lambda(i)*V(:, i);
 
end

