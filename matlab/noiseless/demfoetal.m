%clear all;
HOME = getenv('HOME');
addpath([HOME '/mlprojects/gca/matlab'])
rand('seed', 3.15e5)
randn('seed', 3.15e5);
display = 2;
load 'c:\datasets\ICA\ECG Data\foetal_ecg.dat'



X = foetal_ecg(:, 2:9);
ndata = size(X, 1);

% X = zero meaned data;
meanX = mean(X, 1);
for i = 1:ndata
  X(i, :)  = X(i, :) - meanX;
end

NU_GAUSS = 100;
LATENTDIM = 8;
DATADIM = size(X, 2);
NDATA = size(X, 1);
min_tau = 2.5;%0001;

W = randn(LATENTDIM, DATADIM)*0.1;
[U, V] = eig(cov(X));
NU_TAU = 5*ones(1, LATENTDIM);
SIGMA2_TAU = (NU_TAU-2)./NU_TAU;
W = U*sqrt(diag(1./diag(V)));
W = W';
w = W(:)';
w = [w sqrt(NU_TAU - min_tau)];
options = foptions;
options(1) = 1;
options(9) = 1;
options(14) = 1000;

w = quasinew('objective', w, options, 'gradient', X, min_tau);

W = reshape(w(1:LATENTDIM*LATENTDIM), LATENTDIM, size(X, 2));
SBAR = X*W';

A = inv(W);
NU_TAU = min_tau + w(LATENTDIM*LATENTDIM+1:end).^2;
SIGMA2_TAU = (NU_TAU-2)./NU_TAU;

PCS = find(NU_TAU>NU_GAUSS);

orthoA = A;
orthoA(:, PCS) = orthogonalise(A(:, PCS)')'

orthoW = inv(orthoA);
orthoS = X*orthoW';


