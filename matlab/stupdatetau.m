function [nubar_tau, sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu_tau, sigma2_tau);

if nargin < 4
  sigma2_tau = 1;
end

if nargin < 3
  nu_tau = .5;
end

a = nu_tau/2;
b = a.*sigma2_tau;

latentDim = size(sbar, 2);
ndata = size(sbar, 1);

if length(a) == 1
  a = nrepmat(a, 2, latentDim);
end
if length(b) == 1
  b = nrepmat(b, 2, latentDim);
end

bbar = zeros(ndata, 1);
abar = 0.5+a;

for j = 1:latentDim
%  cbar(j) = 0.5+c(j);
%  dbar(:, j) = (d(j)+0.5*(sbar(:, j).*sbar(:, j) + squeeze(Sigma_s(j, j, :))));
  bbar = b(j) +0.5*(sbar(:, j).*sbar(:, j) + squeeze(Sigma_s(j, j, :)));
  sigma2bar_tau(:, j) = bbar./abar(j);
end

nubar_tau = abar*2;




