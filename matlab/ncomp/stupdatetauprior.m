function stupdatetauprior(method, min_tau)

global LATENTDIM
global NDATA

global NU_TAU
global SIGMA2_TAU

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global SBAR
global SIGMA_S

if nargin < 2
  min_tau = 2.1;
end
if nargin < 1
  method = 'scg';
end
options = foptions;
options(1) = 0;
options(9) = 0;
options(14) = 1000;

switch method
  
 case 'scg' 
   for k = 1:LATENTDIM
     gamma = sqrt(NU_TAU(k)-min_tau);
     param = gamma;
     gamma = scg('stpriorobjective', ...
		    gamma,...
		    options, ...
		    'stpriorgradient', ...
		    TAU, ...
		    LNTAU, ...
		    min_tau, k);
     if gamma*gamma < eps
       gamma = sqrt(eps);
     end
     NU_TAU(k) = min_tau+gamma*gamma;
     SIGMA2_TAU(k) = (NU_TAU(k) - 2)/NU_TAU(k);
   end
 
 case 'quasinew' 
   for k = 1:LATENTDIM
     gamma = sqrt(NU_TAU(k)-min_tau);
     param = gamma;
     gamma = quasinew('stpriorobjective', ...
		    gamma,...
		    options, ...
		    'stpriorgradient', ...
		    TAU, ...
		    LNTAU, ...
		    min_tau, k);
     if gamma*gamma < eps
       gamma = sqrt(eps);
     end
     NU_TAU(k) = min_tau+gamma*gamma;
     SIGMA2_TAU(k) = NU_TAU(k)/(NU_TAU(k) - 2);
   end
end








