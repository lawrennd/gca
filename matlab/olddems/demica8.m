rand('seed', 1e5)
randn('seed', 1e5);
ndata = 400;
niters = 30;

trueLatentDim = 4;
expectedLatentDim = 6
dataDim = 6;

trueA = [round(randn(dataDim, trueLatentDim)*40)/10];
trueBeta = 100;

tauSample = stgamrnd(1, 1, ndata, trueLatentDim);
sSample = randn(ndata, trueLatentDim)./sqrt(tauSample);
X = sSample*trueA' + randn(ndata, dataDim)/sqrt(trueBeta);
%pre whiten data
[q, s] = eig(cov(X));
s = diag(s);
s(find(s<0)) = eps;
U = q*diag(1./sqrt(s));
Uinv = diag(sqrt(s))*q';
s = diag(s);
XU = X*U;% + randn(ndata, dataDim)/sqrt(trueBeta);

trueMix = (trueA'*U)';

%X = X./nrepmat(std(X), 1, size(X, 1));

plot(XU(:, 1), XU(:, 2), 'rx')

% X = zero meaned data;

ndata = size(XU, 1);
dataDim = size(XU, 2);
%latentDim = dataDim;

%V = mixing Weights - Matrix of latentDim Rows, dataDim Columns;
if dataDim > expectedLatentDim
  V = eye(dataDim);
else
  V = eye(expectedLatentDim);
end
V = V(1:dataDim, 1:expectedLatentDim);

% tau = latent noise precision - Matrix of ndata Rows, expectedLatentDim Columns;
tau = ones(ndata, expectedLatentDim);
tol = 1e-3;
% beta = output Noise - Scalar;
beta = 1;%ones(dataDim , 1)*1;
counter = 1;
Vstore = V(:)';
expD = eye(expectedLatentDim);
Dstore = diag(expD)';
betaStore = beta(:)';
alpha = ones(1, expectedLatentDim);
expD2 = expD.*expD;
Sigma = U'*beta*U;
%Sigma = eye(dataDim)*beta;
%Sigma =diag(beta);
for i = 1:niters
  maxchange = tol + 1;
  counter1 = 0;
  while maxchange > tol & counter1 < 100
    oldW = expD*V';
    [sbar, Sigma_s] = updatelatent(XU, expD*V', tau, Sigma);
    %V = updateV2(sbar, Sigma_s, XU, expD, V);
    %[V, expD] = updateV3(sbar, Sigma_s, XU, expD, beta);
    V = updateV(sbar, Sigma_s, XU);
    expD = updateD(sbar, Sigma_s, V, XU);
    expD2 = expD.*expD;
    %[expD expD2] = updateD2(sbar, Sigma_s, V, XU, beta, alpha);
    maxchange = max(max(abs(expD*V' -oldW)));
    fprintf(1, 'Iteration %i max change: %d\n', counter1, maxchange)
    %disp(V'*V)
    counter1 = counter1 + 1;
  end
  %beta = updatebeta2(sbar, Sigma_s, XU, expD, expD2, V); 
  %Sigma = eye(dataDim)*beta;
  beta = updatebeta4(sbar, Sigma_s, XU, expD, expD2, V, Uinv, U); 
  Sigma = U'*beta*U;
  %disp(beta)

  counter = counter + 1;
  Vstore(counter, :) = V(:)';
  Dstore(counter, :) = diag(expD)';
  betaStore(counter, :) = beta(:)';
  tau = stupdatetau(sbar, Sigma_s);
  alpha = updatealpha(expD2);
  %tau = mackay_updatetau(sbar, Sigma_s, tau);
  %disp(tau)
  subplot(3, 1, 1)
  plot(Vstore)
  subplot(3, 1, 2)
  plot(Dstore)
  subplot(3, 1, 3)
  plot(log10(betaStore))
  drawnow
end





W = V'*3;
figure
for i = 1:dataDim
  for j = i+1:dataDim
    subplot(dataDim, dataDim, (i-1)*(dataDim)+ j)
    plot(XU(:, i), XU(:, j), 'rx')
    line([ 0 W(i, i)], [0 W(j, i)])
    line([ 0 W(i, j)], [0 W(j, j)])
    axis image
    xlabel(['Dimension ' num2str(i)])
    ylabel(['Dimension ' num2str(j)])
  end
end




