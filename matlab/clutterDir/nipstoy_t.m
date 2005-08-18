% Toy problem for NIPS paper.


ndata = 1000;


truenu = [2 2 2 100 100 100];
dataDim = length(truenu);
latentDim = 6;
trueA = [round(randn(dataDim, latentDim)*40)/10];
trueBeta = 10000;

% Choose a projection
tempMat = randn(dataDim);
[U, V] = eig(tempMat*tempMat');
project = U;

% Sample a Student-t data set

tX = zeros(ndata, dataDim);
% Sample the taus

tauSample = zeros(ndata, dataDim);	
for i = 1:dataDim
  tauSample(:, i) = stgamrnd(truenu(i), 1, ndata, 1);
end

tX = randn(ndata, dataDim)./sqrt(tauSample);
tX = tX*trueA' + randn(ndata, dataDim)/sqrt(trueBeta);

[At, betat, nut] = learnica(tX, 6);

save nips_gca_t_results
subplot(1, 2, 1), matplot(trueA(:, [1 2 3]))
axis off
axis equal
subplot(1, 2, 2), matplot(At(:, [1 4 6]))
axis off
axis equal
%print -deps z:\tex\projects\ica\diagrams\icarecovered.eps
