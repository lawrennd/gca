function updateA

DATADIM = size(X, 2);
LATENTDIM = size(SBAR, 2);
NDATA = size(X, 1);


sum_ssT = zeros(LATENTDIM);
for n = 1:NDATA
  sum_ssT = sum_ssT + SBAR(n, :)'*SBAR(n, :) + SIGMA_S(:, :, n);
end
U = chol(sum_ssT);
Uinv = eye(LATENTDIM)/U;
inv_sum_ssT = Uinv*Uinv'; 

for i = 1:DATADIM
  A(i, :) = X(:, i)'*SBAR*inv_sum_ssT;
end
%A = A./nrepmat(sum(A.*A, 1), 1, LATENTDIM);