function updatebeta


global NDATA
global DATADIM

global X

global BETA

global A
global SBAR
global SIGMA_S

expectedOutput = (X - SBAR*A');

BETA = 0;
for i = 1:DATADIM
  BETA = BETA + expectedOutput(:, i)'*expectedOutput(:, i);
  for n = 1:NDATA
    BETA = BETA + A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
  end
end
BETA = NDATA*DATADIM/BETA;
