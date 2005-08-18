function [expD, expD2] = updateD2(sbar, Sigma_s, V, UX, beta, alpha)

dataDim = size(UX, 2);
latentDim = size(sbar, 2);
ndata = size(UX, 1);
D = zeros(dataDim);
sigma_D = zeros(dataDim);
sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_s(:, :, n);
end

for j = 1:dataDim
  sigma_D(j, j) = 1/(beta*V(:, j)'*sum_ssT*V(:, j) + alpha(j));
  expD(j, j) = sigma_D(j, j)*beta*UX(:, j)'*sbar*V(:, j);
end
expD2 = expD.*expD + sigma_D;