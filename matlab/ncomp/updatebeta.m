function updatebeta


global NDATA
global DATADIM

global XULINV
global LAMBDA

global BETA

global V
global SBAR
global SIGMA_S

expectedOutput = (XULINV - SBAR*V);

BETA = 0;
for i = 1:DATADIM
  BETA = BETA + expectedOutput(:, i)'*expectedOutput(:, i)*(LAMBDA(i)*LAMBDA(i));
  for n = 1:NDATA
    BETA = BETA + V(:, i)'*SIGMA_S(:, :, n)*V(:, i)*(LAMBDA(i)*LAMBDA(i));
  end
end
BETA = NDATA*DATADIM/BETA;
