function updatebeta


global NDATA
global DATADIM

global X

global BETA

global A
global SBAR
global SIGMA_S

global FANOISE
expectedOutput = (X - SBAR*A');

if FANOISE % Factor analysis type noise model
  BETA = zeros(1, DATADIM);
  for i = 1:DATADIM
    BETA(i) = BETA(i) + expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:NDATA
      BETA(i) = BETA(i) + A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
    end
    BETA(i) = NDATA/BETA(i);
  end
else
  BETA = 0;
  
  for i = 1:DATADIM
    BETA = BETA + expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:NDATA
      BETA = BETA + A(i, :)*SIGMA_S(:, :, n)*A(i, :)';
    end
  end
  BETA = NDATA*DATADIM/BETA;
end