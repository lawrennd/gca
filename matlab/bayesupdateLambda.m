function [lambda, sigma2_lambda] = bayesupdateLambda(sbar, Sigma_s, XU, V, beta, alpha)


dataDim = size(XU, 2);
latentDim = size(V, 2);
ndata = size(XU, 1);
lambda = zeros(1, dataDim);
sum_ssT = sbar'*sbar + sum(Sigma_s, 3);
for i = 1:dataDim
  invsigma2_lambda = alpha(i);
  invsigma2_lambda = invsigma2_lambda + beta(i)*V(:, i)'*sum_ssT*V(:, i);
  sigma2_lambda(i) = 1/invsigma2_lambda;
  for n = 1:ndata
    lambda(i) = lambda(i) + V(:, i)'*sbar(n, :)'*XU(n, i);
  end
  lambda(i) = sigma2_lambda(i)*beta(i)*lambda(i);
end
