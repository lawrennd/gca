function [sbar] = updatesbar(Sigma_s, XUL, V, tau, beta, lambda, choose) 

ndata = size(XUL, 1);
if nargin < 7
  choose = 1:ndata;
end
latentDim = size(tau, 2);
dataDim = size(V, 1);
sbar = zeros(ndata, latentDim);
B = diag(beta*lambda.*lambda);
VBVT = V*B*V';
for n = choose
  sbar(n, :) = XUL(n, :)*B*V'*Sigma_s(:, :, n);
end
