rand('seed', 1e5)
randn('seed', 1e5);
ndata = 10000;
niters = 100;

trueLatentDim = 2;
expectedLatentDim = 2
dataDim = 2;

trueA = [round(randn(dataDim, trueLatentDim)*40)/10];
trueBeta = 10000;

tauSample = stgamrnd(1.5, 1, ndata, trueLatentDim);
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
  V = randn(dataDim);
else
  V = randn(expectedLatentDim);
end
V = V(1:dataDim, 1:expectedLatentDim);
for j = 1:expectedLatentDim
  V(j, :) = V(j, :)/sqrt(V(j, :)*V(j, :)');
end

% tau = latent noise precision - Matrix of ndata Rows, expectedLatentDim Columns;
nu = ones(1, expectedLatentDim);
sigma2_tau = ones(1, expectedLatentDim);
tau = ones(ndata, expectedLatentDim);
tol = 1e-3;
% beta = output Noise - Scalar;
beta = ones(dataDim , 1)*1000;
counter = 1;
Vstore = V(:)';
betaStore = beta(:)';
alpha = ones(1, expectedLatentDim);
Sigma =diag(beta);

for i = 1:niters
  maxchange = tol + 1;
  counter1 = 0;
  while maxchange > tol & counter1 < 100
    oldV = V';
    [sbar, Sigma_s] = updatelatent(XU, V', tau, Sigma);
    V = updateV(sbar, Sigma_s, XU);
    maxchange = max(max(abs(V' -oldV)));
    fprintf(1, 'Iteration %i max change: %d\n', counter1, maxchange)
    counter1 = counter1 + 1;
  end
  %disp(V'*V)
  beta = updatebeta3(sbar, Sigma_s, XU, V); 
  Sigma = diag(beta);

  counter = counter + 1;
  Vstore(counter, :) = V(:)';
  betaStore(counter, :) = beta(:)';
  if ~rem(i, 20)
    [tau, lntau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    [nu, sigma2_tau] = update_tauprior(nu, sigma2_tau, tau, ...
					   lntau, 'fullscg');
  else
    tau = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
  end
  
  %tau = mackay_updatetau(sbar, Sigma_s, tau);
  subplot(3, 1, 1)
  plot(Vstore)
  subplot(3, 1, 2)
  subplot(3, 1, 3)
  plot(log10(betaStore))
  drawnow
end





W = Uinv'*V'*3;
figure
for i = 1:dataDim
  for j = i+1:dataDim
    subplot(dataDim, dataDim, (i-1)*(dataDim)+ j)
    plot(X(:, i), X(:, j), 'rx')
    line([ 0 W(i, i)], [0 W(j, i)])
    line([ 0 W(i, j)], [0 W(j, j)])
    axis image
    xlabel(['Dimension ' num2str(i)])
    ylabel(['Dimension ' num2str(j)])
  end
end
figure
for i = 1:expectedLatentDim
  for j = i+1:expectedLatentDim
    subplot(expectedLatentDim, expectedLatentDim, (i-1)*(expectedLatentDim)+ j)
    plot(sbar(:, i), sbar(:, j), 'rx')
    axis image
    xlabel(['Dimension ' num2str(i)])
    ylabel(['Dimension ' num2str(j)])
  end
end




