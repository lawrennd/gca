function updateV2

global NDATA
global DATADIM
global LATENTDIM

global XULINV

global V
global SBAR
global SIGMA_S

global NU_TAU
global SIGMA2_TAU

sum_ssT = zeros(LATENTDIM);
for n = 1:NDATA
  sum_ssT = sum_ssT + SBAR(n, :)'*SBAR(n, :) + SIGMA_S(:, :, n);
end
C = chol(sum_ssT);
Cinv = eye(LATENTDIM)/C;
inv_sum_ssT = Cinv*Cinv'; 


V = (XULINV'*SBAR*inv_sum_ssT)';


for j = 1:LATENTDIM
  V(j, :) = V(j, :)/sqrt(V(j, :)*V(j, :)');
end

for j = 1:LATENTDIM
  if NU_TAU(j) > 2.1
    V(j, :) = V(j, :)*(NU_TAU(j)-2)/(NU_TAU(j)*SIGMA2_TAU(j));
  else
    V(j, :) = V(j, :)*0.1/(2.1*SIGMA2_TAU(j));

  end
end