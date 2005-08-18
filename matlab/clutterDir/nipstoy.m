% Toy problem for NIPS paper.


ndata = 1000;


sdVector = 10:-1:1;
dataDim = length(sdVector);

% Choose a projection
tempMat = randn(dataDim);
[U, V] = eig(tempMat*tempMat');
project = U;

% Sample a Gaussian data set
gaussX = zeros(ndata, dataDim);
for i = 1:dataDim;
  gaussX(:, i) = randn(ndata, 1)*sqrt(sdVector(i));
end
gaussX = gaussX*project';


[Agauss, betaGauss, nuGauss] = learnica(gaussX, 5);


% Sample a Student-t data set
nu = 2.5;
sigma2 = ((sdVector)*nu-2)/nu;

tX = zeros(ndata, dataDim);
% Sample the taus

tauSample = zeros(ndata, dataDim);	
for i = 1:dataDim
  tauSample(:, i) = stgamrnd(nu, sigma2(i), ndata, 1);
end

tX = randn(ndata, trueLatentDim)./sqrt(tauSample);
tX = tX*project';

[At, betat, nut] = learnica(tX, 5);

