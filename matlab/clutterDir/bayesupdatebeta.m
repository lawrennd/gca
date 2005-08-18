function [abar, bbar] = bayesupdatebeta(sbar, Sigma_s, XULinv, V, Sigma_V, a, b, lambda);

ndata = size(XULinv, 1);
dataDim = size(XULinv, 2);
beta = zeros(1, dataDim);

if nargin < 7
  lambda = ones(1, dataDim);
end
if nargin < 6
  b = 1e-3;
end
if nargin < 5
  a = 1e-3;
end


expectedOutput = (XULinv - sbar*V);
beta = zeros(1, dataDim);
abar = a + 0.5*ndata*dataDim;

bbar = 2*b;
for i = 1:dataDim
  bbar = bbar + expectedOutput(:, i)'*expectedOutput(:, i)*(lambda(i)*lambda(i)) ...
	 ;
  for n = 1:ndata
    bbar = bbar + (V(:, i)'*Sigma_s(:, :, n)*V(:, i) ...
	   + trace(Sigma_s(:, :, n)*Sigma_V(:, :, i)) ...
	   + sbar(n, :)*Sigma_V(:, :, i)*sbar(n, :)')*(lambda(i)*lambda(i));
  end
end
bbar = bbar/2;

