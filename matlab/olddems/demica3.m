rand('seed', 1e5)
randn('seed', 1e5);
ndata = 600;
niters = 100;

latentDim = 4;
dataDim = 4;

trueA = [round(randn(dataDim, latentDim)*40)/10];
trueBeta = 10;

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

%V = mixing Weights - Matrix of latentDim Rows, dataDim Columns;
if dataDim > latentDim
  V = eye(dataDim);
else
  V = eye(latentDim);
end
V = V(1:dataDim, 1:latentDim);

% tau = latent noise precision - Matrix of ndata Rows, latentDim Columns;
tau = ones(ndata, latentDim);
tol = 1e-3;
% beta = output Noise - Scalar;
beta = 10;
counter = 1;
Vstore = V(:)';
betaStore = beta;
alpha = ones(1, latentDim)*0.01;
expD = eye(latentDim);
expD2 = expD.*expD;
for i = 1:niters
  maxchange = tol + 1;
  counter1 = 0;
  while maxchange > tol 
    oldV = V;
    [sbar, Sigma_s] = updatelatent(whiteX, expD*V', tau, beta);
    %      V = updateV2(sbar, Sigma_s, whiteX, D, V);
    %      V = updateV3(sbar, Sigma_s, whiteX, D);
    V = updateV(sbar, Sigma_s, whiteX);
    expD = updateD(sbar, Sigma_s, V, whiteX);
    maxchange = max(max(abs(V -oldV)));
    fprintf(1, 'Iteration %i max change: %d\n', counter1, maxchange)
    %disp(V'*V)
    counter1 = counter1 + 1;
  end
  beta = updatebeta(sbar, Sigma_s, whiteX, expD*V'); 
  %disp(beta)

  counter = counter + 1;
  Vstore(counter, :) = V(:)';
  betaStore(counter) = beta;
  tau = updatetau(sbar, Sigma_s);
  %tau = mackay_updatetau(sbar, Sigma_s, tau);
  %disp(tau)
  subplot(2, 1, 1)
  plot(Vstore)
  subplot(2, 1, 2)
  plot(log10(betaStore))
  drawnow
end





W = V'*3;
figure
for i = 1:dataDim
  for j = i+1:dataDim
    subplot(dataDim, dataDim, (i-1)*(dataDim)+ j)
    plot(whiteX(:, i), whiteX(:, j), 'rx')
    line([ 0 W(i, i)], [0 W(j, i)])
    line([ 0 W(i, j)], [0 W(j, j)])
    axis image
    xlabel(['Dimension ' num2str(i)])
    ylabel(['Dimension ' num2str(j)])
  end
end