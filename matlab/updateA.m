function A = updateA(model, X)

% UPDATEA Update the mixing matrix values.

% GCA

model.dataDim = size(X, 2);
model.latentDim = size(model.sBar, 2);
model.numData = size(X, 1);

A = zeros(model.dataDim, model.latentDim);

sum_ssT = zeros(model.latentDim);
for n = 1:model.numData
  sum_ssT = sum_ssT + model.sBar(n, :)'*model.sBar(n, :) + model.Sigma_s(:, :, n);
end
U = chol(sum_ssT);
Uinv = eye(model.latentDim)/U;
inv_sum_ssT = Uinv*Uinv'; 
for i = 1:model.dataDim
  A(i, :) = X(:, i)'*model.sBar*inv_sum_ssT;
end
