function beta = updatebeta3(sbar, Sigma_s, X, V, Sigma_V);

a = 1e-3;
b = 1e-3;
ndata = size(X, 1);
dataDim = size(X, 2);
beta = zeros(1, dataDim);
% Need to use expected values for A!!!!

for i = 1:dataDim
  for n = 1:ndata
    beta(i) = beta(i) + (X(n, i) - V(:, i)'*sbar(n, :)')^2 + V(:, i)'*Sigma_s(:, :, n)*V(:, i);
    
  end
  beta(i) =  (a + 0.5*ndata)/(b+0.5*beta(i));
end
%beta = (a + 0.5*dataDim*ndata)/(b+0.5*beta);

