function updatelatent

global NDATA
global LATENTDIM

global X
global BETA

global TAU

global A
global SBAR
global SIGMA_S

SBAR = zeros(NDATA, LATENTDIM);
SIGMA_S = zeros(LATENTDIM, LATENTDIM, NDATA);
B = diag(BETA);
ATBA = A'*A*BETA;
for n = 1:NDATA
  invSigma_s = diag(TAU(n, :)) + ATBA;
  C = chol(invSigma_s);
  Cinv = eye(LATENTDIM)/C;
  SIGMA_S(:, :, n) = Cinv*Cinv'; 
  SBAR(n, :) = BETA*X(n, :)*A*SIGMA_S(:, :, n);
end

