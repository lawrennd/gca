function [V, d] = gsupdateV(sbar, Sigma_s, XUL)

dataDim = size(XUL, 2);
latentDim = size(sbar, 2);
ndata = size(XUL, 1);
d = ones(1, latentDim);


sum_ssT = zeros(latentDim);
for n = 1:ndata
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_s(:, :, n);
end
C = chol(sum_ssT);
Cinv = eye(latentDim)/C;
inv_sum_ssT = Cinv*Cinv'; 


for i = 1:dataDim
  V(:, i) = (XUL(:, i)'*sbar*inv_sum_ssT)';
end
[q, d] = eig(V*V');
V = q*diag(sqrt(1./diag(d)))*q'*V;

