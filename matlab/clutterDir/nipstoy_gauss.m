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


Abar = Agauss./nrepmat(sqrt(sum(Agauss.*Agauss)), 1, 10);

gaussX = gaussX - nrepmat(mean(gaussX), 1, 1000);
[U I] = eig(gaussX'*gaussX);
P = Abar'*U;
projections = sqrt(sum(P.^2));

save nips_gca_gaussian_results
