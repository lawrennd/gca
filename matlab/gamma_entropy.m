function h = gamma_entropy(a, b);

% GAMMA_ENTROPY Entropy of a gamma with a gamma prameterised so that E(x) = a/b

% GCA

a = a(:);
b = b(:);

if length(a) ~= 1
  if length(a) ~= length(b)
    error('Gamma parameters must have the same dimensions');
  end
  
  h = sum(log(b)) + sum(digamma(a).*(a-1)) - sum(a) - sum(gammaln(a));
else
  if length(b) == 1
    h = log(b) + digamma(a).*(a-1) - a - gammaln(a);
  else
    h = sum(log(b)) + length(b)*(digamma(a).*(a-1) - a - gammaln(a));
  end
end
h = -h;