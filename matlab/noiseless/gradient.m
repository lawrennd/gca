function lg = gradient(w, X, min_tau);



dataDim = size(X, 2);
latentDim = dataDim;
W = reshape(w(1:latentDim*dataDim), latentDim, dataDim);
gamma = w(latentDim*dataDim + 1:end);
nData = size(X, 1);
V= inv(W);
S = X*W';
W_grad = nData*V';

for j = 1:latentDim
  for i = 1:dataDim
    W_grad(j, i) = W_grad(j, i) + sum(X(:, i).*gradst(S(:, j), gamma(j), min_tau));
  end
end
lg = -W_grad(:)';

gg = zeros(1, latentDim);
for i = 1:latentDim
    gg(i) = gg(i) + 2*gamma(i)*sum(gradgamma(S(:, i), gamma(i), min_tau));
end
gg = -gg;
lg = [lg gg];

function lg = gradst(vals, gamma, min_tau);

nu = min_tau + gamma*gamma;
lg = - (nu+1)*vals/(nu-2).*1./(1+1/(nu-2)*vals.*vals);


function gg = gradgamma(vals, gamma, min_tau);

nu = min_tau + gamma*gamma;
% Mackay's approximation gg = .5*(2/nu + log(nu/(nu+1))) - .5*1/(nu-2) ...
gg = .5*psi((nu+1)/2) - .5*psi(nu/2) - .5*1/(nu-2) ...
      - .5*log(1+vals.*vals/(nu-2)) ...
     + (nu+1)/2*1./(1+vals.*vals/(nu-2)).*vals.*vals/((nu-2)*(nu-2));







