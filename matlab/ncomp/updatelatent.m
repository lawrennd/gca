function updatelatent

global NDATA
global LATENTDIM

global XULINV
global LAMBDA

global BETA

global TAU

global V
global SBAR
global SIGMA_S

SBAR = zeros(NDATA, LATENTDIM);
SIGMA_S = zeros(LATENTDIM, LATENTDIM, NDATA);
B = diag(BETA*LAMBDA.*LAMBDA);
VBVT = V*B*V';
for n = 1:NDATA
  invSigma_s = diag(TAU(n, :)) + VBVT;
  C = chol(invSigma_s);
  Cinv = eye(LATENTDIM)/C;
  SIGMA_S(:, :, n) = Cinv*Cinv'; 
  SBAR(n, :) = XULINV(n, :)*B*V'*SIGMA_S(:, :, n);
end

