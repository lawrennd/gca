rand('seed', 1e5)
randn('seed', 1e5);
ndata = 300;
niters = 100;

latentDim = 4;
dataDim = 4;

trueA = [round(randn(dataDim, latentDim)*40)/10];
trueBeta = 1;

tauSample = gamrnd(0.5, 0.5, ndata, latentDim);
sSample = randn(ndata, latentDim)./sqrt(tauSample);
X = sSample*trueA' + randn(ndata, dataDim)/sqrt(trueBeta);
%pre whiten data
[q, s] = eig(cov(X));
whiteningMatrix = q*diag(1./sqrt(diag(s)));
whiteX = X*whiteningMatrix;

trueMix = (trueA'*whiteningMatrix)';

%X = X./nrepmat(std(X), 1, size(X, 1));

plot(whiteX(:, 1), whiteX(:, 2), 'rx')

% X = zero meaned data;

ndata = size(whiteX, 1);
dataDim = size(whiteX, 2);
%latentDim = dataDim;

%A = mixing Weights - Matrix of dataDim Rows, latentDim Columns;
A = randn(dataDim, latentDim);

% tau = latent noise precision - Matrix of ndata Rows, latentDim Columns;
tau = rand(ndata, latentDim);

% beta = output Noise - Scalar;
beta = 1;
%tau = tauSample;
%beta = trueBeta;
%A = trueA;
%sbar = sSample;
%[sbar, Sigma_beta_s] = updatelatent(X, A, tau, beta);
counter = 1;
Astore = A(:)';
betaStore = beta;
for i = 1:niters
  for j = 1:1
    for k = 1:10
      [sbar, Sigma_beta_s] = updatelatent(whiteX, A, tau, beta);
      A = updateA(sbar, Sigma_beta_s, whiteX);
    
    end
    beta = updatebeta(sbar, Sigma_beta_s, whiteX, A); 
    %disp(beta)
  end
  counter = counter + 1;
  Astore(counter, :) = A(:)';
  betaStore(counter) = beta;
  tau = mackay_updatetau(sbar, Sigma_beta_s, tau);
  %disp(tau)
  subplot(2, 1, 1)
  plot(Astore)
  subplot(2, 1, 2)
  plot(log10(betaStore))
  drawnow
end








