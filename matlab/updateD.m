function D = updateD(sbar, Sigma_s, V, XU)

dataDim = size(XU, 2);
latentDim = size(sbar, 2);
ndata = size(XU, 1);
D = zeros(dataDim);

sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_s(:, :, n);
end
for j = 1:dataDim
  D(j, j) = 1/(V(:, j)'*sum_ssT*V(:, j))*XU(:, j)'*sbar*V(:, j);
end