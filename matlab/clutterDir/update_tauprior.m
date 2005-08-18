function [nu_tau, sigma2_tau] = update_tauprior(nu_tau, sigma2_tau, expTau, expLnTau, method)


if nargin < 5
  method = 'scg';
end
options = foptions;
options(1) = 1;

switch method
 case 'moments'
  sigma2_tau = 1./mean(expTau);
  %  varTau = mean(mean(model.abar_tau./(model.bbar_tau.*model.bbar_tau)));
  %  
  %  model.b_tau = meanTau./varTau;
  %  
  %  model.a_tau = meanTau.*model.b_tau;
  
 case 'scg' 
  lnnu_tau = log(nu_tau);
  lnnu_tau = scg('tauobjective', lnnu_tau, options, 'taugradient',  expTau, ...
		 expLnTau, sigma2_tau);
  nu_tau = exp(lnnu_tau);
 
 case 'quasinew' 
  lnnu_tau = log(nu_tau);
  lnnu_tau = quasinew('tauobjective', lnnu_tau, options, 'taugradient', ...
		      expTau, expLnTau, sigma2_tau);
  nu_tau = exp(lnnu_tau);

 case 'fullscg'
  lnparam = log([nu_tau sigma2_tau]);
  lnparam = scg('tauobjective', lnparam, options, 'taugradient', expTau, ...
		expLnTau);
  latentDim = size(expTau, 2);
  nu_tau = exp(lnparam(1:latentDim));
  sigma2_tau = exp(lnparam(latentDim + 1:end));
 
 case 'fullquasinew'
  lnparam = log([nu_tau sigma2_tau]);
  lnparam = quasinew('tauobjective', lnparam, options, 'taugradient', ...
		     expTau, expLnTau);
  latentDim = size(expTau, 2);
  nu_tau = exp(lnparam(1:latentDim));
  sigma2_tau = exp(lnparam(latentDim + 1:end));
  
end








