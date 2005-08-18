function [V, d] = updateV(sbar, Sigma_s, XUL, choose)


dataDim = size(XUL, 2);
latentDim = size(sbar, 2);
ndata = size(XUL, 1);

if nargin < 4
  choose = 1:ndata;
end
d = ones(1, latentDim);


sum_ssT = zeros(latentDim);
for n = choose
  sum_ssT = sum_ssT + sbar(n, :)'*sbar(n, :) + Sigma_s(:, :, n);
end
C = chol(sum_ssT);
Cinv = eye(latentDim)/C;
inv_sum_ssT = Cinv*Cinv'; 


for i = 1:dataDim
  V(:, i) = (XUL(choose, i)'*sbar(choose, :)*inv_sum_ssT)';
end
