function beta = updateB(sbar, Sigma_s, XULinv, V, lambda, U);

ndata = size(XULinv, 1);
dataDim = size(XULinv, 2);

expectedOutput = (XULinv - sbar*V)*diag(lambda)*U';

beta = 0;
for i = 1:dataDim
  beta = beta + expectedOutput(:, i)'*expectedOutput(:, i)*(lambda(i)*lambda(i));
  for n = 1:ndata
    beta = beta + U(lambda(i)*V(:, i)'*Sigma_s(:, :, n)*V(:, i)*lambda(i));
  end
end
beta = ndata*dataDim/beta;
function B = updateB(sbar, Sigma_s, X, V);
A = V';
a = 0.5;
b = 0.5;
ndata = size(X, 1);
dataDim = size(X, 2);
beta = zeros(1, dataDim);
% Need to use expected values for A!!!!
output = (X - sbar*A');
B = zeros(dataDim);
for n = 1:ndata
  B = B + output(n, :)'*output(n, :) + A*Sigma_s(:, :, n)*A';
end
B = inv(B/ndata);
