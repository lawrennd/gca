function stupdatetau

global NDATA
global LATENTDIM

global NU_TAU
global SIGMA2_TAU

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global SBAR
global SIGMA_S

a = NU_TAU/2;
b = a.*SIGMA2_TAU;

if length(a) == 1
  a = nrepmat(a, 2, LATENTDIM);
end
if length(b) == 1
  b = nrepmat(b, 2, LATENTDIM);
end

bbar = zeros(NDATA, 1);
abar = 0.5+a;

for j = 1:LATENTDIM
  bbar = b(j) +0.5*(SBAR(:, j).*SBAR(:, j) + squeeze(SIGMA_S(j, j, :)));
  SIGMA2BAR_TAU(:, j) = bbar./abar(j);
end

NUBAR_TAU = abar*2;




