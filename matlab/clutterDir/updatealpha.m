function alpha = updatealpha(U, V, lambda, sigma2_lambda);

a = 1e-6;
b = 1e-6;

latentDim = size(V, 2);
dataDim = size(V, 1);
alpha = zeros(1, latentDim);
expA = U*diag(lambda)*V';
sigma_lambda = sqrt(sigma2_lambda);
expATA = expA'*expA + U*diag(sigma_lambda)*V'*V*diag(sigma_lambda)*U';
alpha = (0.5*dataDim+a)./(b+0.5*diag(expATA));


  




