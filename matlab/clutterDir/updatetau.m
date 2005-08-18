function [cbar, dbar] = updatetau(sbar, Sigma_s, c, d);

latentDim = size(sbar, 2);
ndata = size(sbar, 1);

tau = zeros(ndata, latentDim);

for j = 1:latentDim
  cbar(j) = 0.5+c(j);
  dbar(:, j) = (d(j)+0.5*(sbar(:, j).*sbar(:, j) + squeeze(Sigma_s(j, j, :))));
end


  




