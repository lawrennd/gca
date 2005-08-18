function [sigma2_tau, nu_tau] = stupdatetauprior(model, X, method, min_tau)

% STUPDATETAUPRIOR Update the values of nu_tau.

% GCA

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

nu_tau = model.nu_tau;
sigma2_tau = model.sigma2_tau;

switch method
  
 case 'scg' 
   for k = 1:model.latentDim
     gamma = sqrt(nu_tau(k)-min_tau);
     param = gamma;
     gamma = scg('stpriorobjective', ...
		    gamma,...
		    options, ...
		    'stpriorgradient', ...
		    model.tau, ...
		     model.lntau, ...
		    min_tau, k);
     if gamma*gamma < eps
       gamma = sqrt(eps);
     end
     nu_tau(k) = min_tau+gamma*gamma;
     sigma2_tau(k) = (nu_tau(k) - 2)/nu_tau(k);
   end
 
 case 'quasinew' 
   for k = 1:model.latentDim
     gamma = sqrt(nu_tau(k)-min_tau);
     param = gamma;
     gamma = quasinew('stpriorobjective', ...
		    gamma,...
		    options, ...
		    'stpriorgradient', ...
		    model.tau, ...
		     model.lntau, ...
		    min_tau, k);
     if gamma*gamma < eps
       gamma = sqrt(eps);
     end
     nu_tau(k) = min_tau+gamma*gamma;
     sigma2_tau(k) = nu_tau(k)/(nu_tau(k) - 2);
   end
end








