function A = updateA(sbar, Sigma_beta_s, X)

dataDim = size(X, 2);
latentDim = size(sbar, 2);
ndata = size(X, 1);


sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_beta_s(:, :, n);
end
U = chol(sum_ssT);
Uinv = eye(latentDim)/U;
inv_sum_ssT = Uinv*Uinv'; 

for i = 1:dataDim
  A(i, :) = X(:, i)'*sbar*inv_sum_ssT;
end
%A = A./nrepmat(sum(A.*A, 1), 1, latentDim);