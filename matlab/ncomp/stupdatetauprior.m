function stupdatetauprior(method, min_tau)

global LATENTDIM
global NDATA

global NU_TAU
global SIGMA2_TAU

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global V
global SBAR
global SIGMA_S

global LAGRANGE

if nargin < 2
  min_tau = 2.1;
end
if nargin < 1
  method = 'scg';
end
options = foptions;
%options(1) = 1;
%options(9) = 1;
options(14) = 1000;

switch method
 case 'moments'
  SIGMA2_TAU = 1./mean(TAU);
  varTau = mean(mean(model.abar_tau./(model.bbar_tau.*model.bbar_tau)));
  a_tau = meanTau./varTau;
  b_tau = meanTau.*model.b_tau;
  
 case 'scg' 
  bndOld = sticabound/NDATA;

  for i = 1:10
    lnNU_TAU = log(NU_TAU);
    lnNU_TAU = scg('sttauobjective', lnNU_TAU, options, 'sttaugradient',  TAU, ...
		   LNTAU, SIGMA2_TAU);
    NU_TAU = exp(lnNU_TAU);
    bnd = sticabound/NDATA;
    disp(bndOld-bnd)
    bndOld = bnd;
    NU_TAU(find(NU_TAU<min_tau)) = min_tau;
    bnd = sticabound/NDATA;
    disp(bndOld - bnd)
    bndOld = bnd;
    SIGMA2_TAU = (NU_TAU-2)./NU_TAU;
    bnd = sticabound/NDATA;
    disp(bndOld - bnd)
    bndOld = bnd;
    LAGRANGE = (SIGMA2_TAU-2)./2 .* (sum(TAU, 1) - NDATA./SIGMA2_TAU);

  end
 
 case 'quasinew' 
  lnNU_TAU = log(NU_TAU);
  lnNU_TAU = quasinew('sttauobjective', lnNU_TAU, options, 'sttaugradient', ...
		      TAU, LNTAU, SIGMA2_TAU);
  NU_TAU = exp(lnNU_TAU);
  NU_TAU(find(NU_TAU<min_tau)) = min_tau;
  SIGMA2_TAU = (NU_TAU-2)./NU_TAU;

 case 'scgconstrained' 
  if isempty(LAGRANGE)
    LAGRANGE = ones(1, LATENTDIM);
  end
  param = [real(log([NU_TAU-min_tau SIGMA2_TAU])) LAGRANGE];
  param = scg('sttauconstrainedobjective', param, options, 'sttauconstrainedgradient', TAU, ...
		LNTAU, min_tau);
  NU_TAU = min_tau+exp(param(1:LATENTDIM));
  NU_TAU(find(isinf(NU_TAU))) = 1000;
  NU_TAU(find(NU_TAU>1000)) = 1000;
  SIGMA2_TAU = exp(param(LATENTDIM + 1:2*LATENTDIM));
  SIGMA2_TAU(find(isinf(SIGMA2_TAU))) = 1000;
  SIGMA2_TAU(find(SIGMA2_TAU>1000)) = 1000;
  SIGMA2_TAU(find(SIGMA2_TAU<eps)) = eps;
  LAGRANGE = param(LATENTDIM*2 + 1:end);

 case 'quasinewconstrained' 
  if isempty(LAGRANGE)
    LAGRANGE = zeros(1, LATENTDIM)*1e3;
  end
  param = [real(log([NU_TAU-min_tau SIGMA2_TAU])) LAGRANGE];
  param = quasinew('sttauconstrainedobjective', param, options, ...
	      'sttauconstrainedgradient', TAU, ...
	      LNTAU, min_tau);
  NU_TAU = min_tau+exp(param(1:LATENTDIM));
  SIGMA2_TAU = exp(param(LATENTDIM + 1:2*LATENTDIM));
  
  LAGRANGE = param(LATENTDIM*2 + 1:end);

 case 'fullscg'
  param = log([NU_TAU SIGMA2_TAU]);
  param = scg('sttauobjective', param, options, 'sttaugradient', TAU, ...
		LNTAU);
  NU_TAU = exp(param(1:LATENTDIM));
  SIGMA2_TAU = exp(param(LATENTDIM + 1:end));
 
 case 'fullquasinew'
  param = log([NU_TAU SIGMA2_TAU]);
  param = quasinew('sttauobjective', param, options, 'sttaugradient', ...
		     TAU, LNTAU);
  NU_TAU = exp(param(1:LATENTDIM));
  SIGMA2_TAU = exp(param(LATENTDIM + 1:end));
  
end








