function stupdatetauprior2(method, min_tau)

global LATENTDIM
global NDATA

global NU_TAU
global SIGMA2_TAU
global LAGRANGE

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global V
global SBAR
global SIGMA_S

if nargin < 2
  min_tau = 2.1;
end

for i = 1:10
  for j = 1:LATENTDIM
    LAGRANGE(j) = (SIGMA2_TAU(j) -2)/2 * (sum(TAU(:, j)) - NDATA/ ...
				     SIGMA2_TAU(j));
  end
  options = foptions;
  options(1) = 1;
  for j = 1:LATENTDIM
    [NU_TAU(j), fval] = fminbnd('nulagrangian', min_tau, 100, [], min_tau, j);
  end
  for j = 1:LATENTDIM
    SIGMA2_TAU(j) = (NU_TAU(j)-2)/NU_TAU(j);
  end
end