function [c, d, obj] = updatetauprior(c, d, expTau, expLnTau, method)


if nargin < 5
  method = 'scg';
end
options = foptions;
%options(1) = 1;
%options(9) = 1;

switch method
 case 'moments'
  meanTau = mean(mean(model.abar_tau./model.bbar_tau));
  varTau = mean(mean(model.abar_tau./(model.bbar_tau.*model.bbar_tau)));
  
  d = meanTau./varTau;
  
  c = meanTau.*d;

 
 case 'scg'
  lnparam = log([c d]);
  lnparam = scg('tauobjective', lnparam, options, 'taugradient', expTau, ...
		expLnTau);
  latentDim = size(expTau, 2);
  c = exp(lnparam(1:latentDim));
  d = exp(lnparam(latentDim + 1:end));
 
 case 'quasinew'
  lnparam = log([c d]);
  lnparam = quasinew('tauobjective', lnparam, options, 'taugradient', expTau, ...
		expLnTau);
  latentDim = size(expTau, 2);
  c = exp(lnparam(1:latentDim));
  d = exp(lnparam(latentDim + 1:end));

end








