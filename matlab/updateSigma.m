function Sigma = updateSigma(sbar, Sigma_s, X, V);
A = V';
a = 0.5;
b = 0.5;
ndata = size(X, 1);
dataDim = size(X, 2);
beta = zeros(1, dataDim);
% Need to use expected values for A!!!!
meanOut =sum(X - sbar*A', 1)';
Sigma = meanOut*meanOut';

