function tau = mackay_updatetau(sbar, Sigma_beta_s, tau, c, d);

latentDim = size(sbar, 2);
ndata = size(sbar, 1);


%tau = zeros(ndata, latentDim);

for j = 1:latentDim
  for n = 1:ndata
    gamma = 1 - tau(n, j)*Sigma_beta_s(j, j, n);
    tau(n, j) = (gamma + 2*c(j))/(sbar(n, j)*sbar(n, j) + 2*d(j));
  end
end


  




