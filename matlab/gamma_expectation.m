function exp = gamma_expectation(a, b, expX, expLnX);

% GAMMA_EXPECTATION Compute the expectation of a Gamma distribution.

% GCA 

expX = expX(:);
expLnX = expLnX(:);

if length(expX) ~= length(expLnX)
  error('Mean and geometric mean must have the same dimensions');
end

a = a(:);
b = b(:);

if length(a) ~= length(b)
  error('Gamma parameters must have the same dimensions');
end

dim = length(expX);

if length(a) == 1
  constPart = dim*(xlogy(a, b) - gammaln(a));
  varPart = (a-1)*sum(expLnX) - b*sum(expX);
else
  if length(a) ~= dim
    error('Gamma parameters must have dimension 1 or the same as expectations')
  else
    constPart = sum(xlogy(a, b)) - sum(gammaln(a));
    varPart = sum((a-1).*expLnX) - sum(b.*expX);
  end
end

exp = constPart + varPart;
