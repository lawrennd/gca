function updateV

global NDATA
global DATADIM
global LATENTDIM

global XULINV

global V
global SBAR
global SIGMA_S

sum_ssT = zeros(LATENTDIM);
for n = 1:NDATA
  sum_ssT = sum_ssT + SBAR(n, :)'*SBAR(n, :) + SIGMA_S(:, :, n);
end
C = chol(sum_ssT);
Cinv = eye(LATENTDIM)/C;
inv_sum_ssT = Cinv*Cinv'; 


for i = 1:DATADIM
  V(:, i) = (XULINV(:, i)'*SBAR(:, :)*inv_sum_ssT)';
end
