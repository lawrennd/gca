function beta = updatebeta(sbar, Sigma_s, XULinv, V, lambda);

ndata = size(XULinv, 1);
dataDim = size(XULinv, 2);

expectedOutput = (XULinv - sbar*V);

beta = 0;
for i = 1:dataDim
  beta = beta + expectedOutput(:, i)'*expectedOutput(:, i)*(lambda(i)*lambda(i));
  for n = 1:ndata
    beta = beta + V(:, i)'*Sigma_s(:, :, n)*V(:, i)*(lambda(i)*lambda(i));
  end
end
beta = ndata*dataDim/beta;
