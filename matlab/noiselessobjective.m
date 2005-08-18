function ll = noiselessobjective(w, X, min_tau)

% NOISELESSOBJECTIVE the objective function for the noiseless GCmodel.A algorithm.

% GCA

dataDim = size(X, 2);
latentDim = dataDim;
nData = size(X, 1);
W = reshape(w(1:latentDim*dataDim), latentDim, dataDim);
gamma = w(latentDim*dataDim + 1: end);

ll = nData*log(det(W));

S = X*W';
for i = 1:latentDim 
  ll = ll + sum(log_studentt(S(:, i), gamma(i), min_tau));
end

ll = -ll;

function lp = log_studentt(vals, gamma, min_tau);

nu = gamma*gamma + min_tau;
sigma2 = (nu-2)/nu;
lp = gammaln((nu+1)/2) - gammaln(nu/2) ...
     - .5*log(nu*pi*sigma2) - (nu+1)/2*log(1+vals.*vals/(nu*sigma2));





