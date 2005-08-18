function [abar_alpha, bbar_alpha] = bayesupdatealpha(V, Sigma_V, a_alpha, ...
						  b_alpha);

dataDim = size(V, 2);
latentDim = size(V, 1);
abar_alpha = a_alpha+dataDim*latentDim/2;
abar_alpha = repmat(abar_alpha, 1, latentDim);
bbar_alpha = zeros(1, latentDim);
for i = 1:dataDim
  bbar_alpha = bbar_alpha + (V(:, i).*V(:, i))' + diag(Sigma_V(:, :, i))';
end
bbar_alpha = sum(bbar_alpha)*0.5 + b_alpha;
bbar_alpha = bbar_alpha*ones(1, latentDim);
