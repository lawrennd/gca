function V = updateV2(sbar, Sigma_beta_s, UX, D, V)

dataDim = size(UX, 2);
latentDim = size(sbar, 2);
ndata = size(UX, 1);
D = diag(D);

sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_beta_s(:, :, n);
end
C = chol(sum_ssT);
Cinv = eye(latentDim)/C;
inv_sum_ssT = Cinv*Cinv'; 

for j = 1:dataDim
  lambda(j) = V(:, j)'*sum_ssT*V(:, j) - 1/D(j) *UX(:, j)'*sbar*V(:, j);
end
V(:, 1) = 1/D(1)*(UX(:, 1)'*sbar*inv(sum_ssT - lambda(1)*eye(latentDim)))';

for j = 2:dataDim
  V(:, j) = 1/D(j)*(UX(:, j)'*sbar*inv(sum_ssT - lambda(j)*inv(eye(latentDim) ...
						  - V(:, 1:j-1)*V(:, 1:j-1)')))';
end
%for j = 1:dataDim
%  
%  V(:, j) = (UX(:, j)'*sbar*inv_sum_ssT)';
%  V(:, j) = V(:, j)/sqrt(V(:, j)'*V(:, j));
%end
%A = A./nrepmat(sum(A.*A, 1), 1, latentDim);