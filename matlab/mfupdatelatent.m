function [sbar, Sigma_s, sbardiff] = mfupdatelatent(XULinv, V, tau, beta, lambda, sbar) 

ndata = size(XULinv, 1);
latentDim = size(tau, 2);
dataDim = size(V, 1);
%sbar = zeros(ndata, latentDim);
Sigma_s = zeros(latentDim, latentDim, ndata);
Bl2 = (beta*lambda.*lambda)';
sbardiff = 0;
V2 = V'.*V';
for j = 1:latentDim
  Bl2V2(:, j) = Bl2.*V2(:, j);
end
sindex = randperm(latentDim);
for n = 1:ndata
  for j = sindex
    notj = [1:j-1 j+1:latentDim];
    oldsbar = sbar(n, j);
    sigma2 = 1/(sum(Bl2V2(:, j))) + tau(n, j);
    sbar(n, j) = sigma2*...
	sum(Bl2.*(XULinv(n, :)' - V(notj, :)'*sbar(n, notj)') ...
	    .*V(j, :)');
    Sigma_s(j, j, n) = sigma2;
    tempsbardiff = abs(sbar(n, j) - oldsbar);
    if tempsbardiff > sbardiff
      sbardiff = tempsbardiff;
    end
  end
end
