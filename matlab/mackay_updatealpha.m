function [alpha, gamma] = mackay_updatealpha(V, Sigma_V);

a = 0.01;
b = 0.01;

latentDim = size(V, 2);
dataDim = size(V, 1);
alpha = zeros(1, latentDim);
sqExp = zeros(1, latentDim);
for i = 1:dataDim
  for j = 1:latentDim
    sqExp(j) = sqExp(j) + Sigma_V(j, j, i);
  end
end
for j = 1:latentDim
  gamma(j) = 1 - alpha(j) * sqExp(j);
  alpha(j) = (gamma(j) + 2*a)/(V(j, :)*V(j, :)' + 2*b);
end




  




