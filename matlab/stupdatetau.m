function [sigma2Bar_tau, nuBar_tau] = stupdatetau(model, X)

% STUPDATETAU Update the values of model.tau.

% GCA

a = model.nu_tau/2;
b = a.*model.sigma2_tau;

if length(a) == 1
  a = nrepmat(a, 2, model.latentDim);
end
if length(b) == 1
  b = nrepmat(b, 2, model.latentDim);
end

bbar = zeros(model.numData, 1);
abar = 0.5+a;

sigma2Bar_tau = zeros(size(model.sigma2Bar_tau));
for j = 1:model.latentDim
  bbar = b(j) +0.5*(model.sBar(:, j).*model.sBar(:, j) + squeeze(model.Sigma_s(j, j, :)));
  sigma2Bar_tau(:, j) = bbar./abar(j);
end

nuBar_tau = abar*2;






