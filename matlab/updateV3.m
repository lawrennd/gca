function [V, beta] = updateV3(sbar, Sigma_s, XU)

dataDim = size(XU, 2);
latentDim = size(sbar, 2);
ndata = size(XU, 1);

sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_s(:, :, n);
end
sum_ssT = sum_ssT;
[V, lambda] =  eig(sum_ssT);
V = V';
beta = 1./diag(lambda);
