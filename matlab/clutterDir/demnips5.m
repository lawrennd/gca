clear all;
rand('seed', 1e3)
randn('seed', 1e3);
ndata = 3500;
niters = 500;
display = 2;
load 'c:\datasets\ICA\ECG Data\foetal_ecg.dat'
expectedLatentDim = 10;

X = foetal_ecg(:, 2:6); %6500)';

%clear MEG_art;

%X = X(randperm(1:1000), :);

ndata = size(X, 1);
dataDim = size(X, 2);

% X = zero meaned data;
meanX = mean(X, 1);
for i = 1:ndata
  X(i, :)  = X(i, :) - meanX;
end


%pre whiten data
[q, s] = eig(cov(X));
s = diag(s);
s(find(s<0)) = eps;
lambda = sqrt(s)';
sigma2_lambda = 1e-6*ones(1, dataDim);
U = q;
UT = q';
s = diag(s);
XULinv = X*U*diag(1./lambda);

plot(XULinv(:, 1), XULinv(:, 2), 'rx')

% X = zero meaned data;

ndata = size(XULinv, 1);
dataDim = size(XULinv, 2);
%latentDim = dataDim;

%V = mixing Weights - Matrix of latentDim Rows, dataDim Columns;
if dataDim > expectedLatentDim
  V = eye(dataDim);
else
  V = eye(expectedLatentDim);
end
V = randn(expectedLatentDim, dataDim);
%for j = 1:expectedLatentDim
%  V(j, :) = V(j, :)/sqrt(V(j, :)*V(j, :)');
%end

% Tolerances
tol = 1e-5;
tautol = 1e-3;
lntautol = 1e-3;
lltol = 1e-3;
 
nu = 1*ones(1, expectedLatentDim);
sigma2_tau = 1*ones(1, expectedLatentDim);

startBeta = 1/mean(diag(cov(X)));
beta = 100*startBeta;



nubar_tau = nrepmat(nu, 1, ndata);
sigma2bar_tau = nrepmat(sigma2_tau, 1, ndata);
tau = 1./sigma2bar_tau;
for j = 1:expectedLatentDim
  lntau(:, j) = psi(nubar_tau(j)/2) ...
      - log(nubar_tau(j)/2) ...
      - log(sigma2bar_tau(:,j));
end

counter = 1;
Vstore = V(:)';
betaStore = beta(:)';
nuStore = nu(:)';
D = eye(expectedLatentDim);
llCounter = 1;
sbar = zeros(ndata, expectedLatentDim);
[sbar, Sigma_s] = updatelatent(XULinv, V, tau, beta, lambda);
lll(llCounter) = sticabound(V, lambda, XULinv, beta, ...
			    nu, sigma2_tau, ...
			    sbar, Sigma_s, ...
			    nubar_tau, sigma2bar_tau);
for iter = 1:niters
  estepCount = 0;
  lntauchange = lntautol + 1;
  tauchange = tautol + 1;
  %  while (tauchange > tautol | lntauchange > lntautol) & estepCount < 200
  %    oldtau = tau;
  %    oldlntau = lntau;
  %    estepCount = estepCount + 1;
  
  maxchange = tol + 1;
  counter1 = 0;
  
  while maxchange > tol & counter1 < 500
    
    oldV = V';
    [sbar, Sigma_s] = updatelatent(XULinv, V, tau, beta, ...
					       lambda);
    V = updateV(sbar, Sigma_s, XULinv);
    [nubar_tau sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    tau = 1./sigma2bar_tau;
    for j = 1:expectedLatentDim
      lntau(:, j) = psi(nubar_tau(j)/2) ...
	  - log(nubar_tau(j)/2) ...
	  - log(sigma2bar_tau(:,j));
    end

    maxchange = max(max(abs(V' -oldV)));
    counter1 = counter1 + 1;
    
    if display > 1
      fprintf(1, 'Iteration %i.%i max V change: %d\n', ...
	      iter, counter1, maxchange)
      drawnow
    end
    if display > 2
      figure(5)
	for i = 1:expectedLatentDim
	  for j = i+1:expectedLatentDim
	    subplot(expectedLatentDim, expectedLatentDim, (i-1)*(expectedLatentDim)+ j)
	    plot(sbar(:, i), sbar(:, j), 'rx')
	    axis image
	    set(gca, 'XTickLabel', '')
	    set(gca, 'YTickLabel', '')
	  end
	end
	drawnow
    end
    %    % End of SV loop
  end
  
  beta = updatebeta(sbar, Sigma_s, XULinv, V, lambda); 
  
  %    tauchange = max(max(abs(tau - oldtau)));
  %    lntauchange = max(max(abs(lntau - oldlntau)));
  %    fprintf(1, 'Iteration %i.%i max tau change: %d, max lntau change %d\n', ...
  %	    iter, estepCount, tauchange, lntauchange)
  counter = counter + 1;
  
  if display
    Vstore(counter, :) = V(:)';
    betaStore(counter, :) = beta(:)';
    nuStore(counter, :) = nu(:)';
  end
  
  if display > 1
    figure(1)
    subplot(4, 1, 1)
    plot(Vstore)
    subplot(4, 1, 2)
    plot(log10(betaStore))
    subplot(4, 1, 3)
    plot(log10(nuStore))
    subplot(4, 1, 4)
    plot(1:llCounter, lll(1:llCounter)/ndata);
    drawnow
  end
  % End of tau loop  
  %  end
  
  
  % M step in latent distributions.

  taucounter = 0;
  entropy_tau = 0;
  for j = 1:expectedLatentDim
    entropy_tau = entropy_tau + stgamma_entropy(nubar_tau(j), sigma2bar_tau(:, j)); 
  end
  obj = sttauobjective(log(nu), tau, lntau, sigma2_tau) - entropy_tau;
  objdif = 1+ 1e-3;
  while objdif > tol & taucounter < 500
    oldobj = obj;
    taucounter = taucounter + 1;
    [nu, sigma2_tau] = stupdatetauprior(nu, sigma2_tau, tau, lntau, 'scg');
    [nubar_tau sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    
    entropy_tau = 0;
    tau = 1./sigma2bar_tau;
    for j = 1:expectedLatentDim
      lntau(:, j) = psi(nubar_tau(j)/2) ...
	  - log(nubar_tau(j)/2) ...
	  - log(sigma2bar_tau(:, j));
      entropy_tau = entropy_tau + stgamma_entropy(nubar_tau(j), sigma2bar_tau(:, j)); 
    end
    
    obj = sttauobjective(log(nu), tau, lntau, sigma2_tau) - entropy_tau;
    objdif = oldobj - obj;
    if display > 1
      if objdif < 0
	warning(['Tau objective increased by ' num2str(objdif)])
      end
      fprintf('Tau Iteration %i.%i, tau change %d\n', iter, taucounter, objdif);
    end
    %disp(objdif);
  end
  llCounter = llCounter + 1;
  lll(llCounter) = sticabound(V, lambda, XULinv, beta, ...
			     nu, sigma2_tau, ...
			     sbar, Sigma_s, ...
			     nubar_tau, sigma2bar_tau);
  
  llldiff = lll(llCounter) - lll(llCounter-1);
  fprintf('Iteration %i, log likelihood change: %d\n', iter, llldiff)
  if llldiff < tol
    break
  end
  A = U*diag(lambda)*V';
  save nips5save A beta nu
end

A = U*diag(lambda)*V';
figure(3)
for i = 1:expectedLatentDim
  for j = i+1:expectedLatentDim
    subplot(expectedLatentDim, expectedLatentDim, (i-1)*(expectedLatentDim)+ j)
    plot(sbar(:, i), sbar(:, j), 'rx')
    axis image
  end
end


figure(4)
clf
yTopStore = 0;
[void, index] = sort(nu);
counter = 0;
SLim(1) = min(min(s));
SLim(2) = max(max(s));
plotHeight = 1/expectedLatentDim-0.05/expectedLatentDim;
for i = index
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (expectedLatentDim - counter)/expectedLatentDim+0.01 0.77 plotHeight]);
  plot(1:size(sbar, 1), sbar(:, i));
  yLim = get(gca, 'YLim');
  if yLim(2) < 4;
    yLim(2) = 4;
  end
    if yLim(1) > -4;
    yLim(1) = -4;
  end
  
  set(gca, 'YLim', [-10 10]);
  set(gca, 'XTickLabel', '')
  set(gca, 'YTick', [yLim(1) 0 yLim(2)])
  axis 
  ax(counter, 2) = axes('position', [0.81 (expectedLatentDim - counter)/expectedLatentDim+0.01 0.08 plotHeight]);
  x = linspace(-4, 4, 200);
  y = tpdf(x, nu(i), sigma2_tau(i));
  plot(x, y, 'b-')
  yLim = get(gca, 'YLim');
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  set(gca, 'XLim', [-10 10]);
  set(gca, 'XTickLabel', '')%axis off
  set(gca, 'YTickLabel', '')
  ax(counter, 3) = axes('position', [0.91 (expectedLatentDim - counter)/expectedLatentDim+0.01 0.08 plotHeight]);
  disp(nu(i))
  index2 = find(sbar(:, i) > -10 & sbar(:, i) < 10);
  index2 = 1:length(index2);
  [heights, centres] = hist(sbar(index2, i), 30);
  heights = heights/sum(heights);
  binwidth = centres(2) - centres(1);
  heights = heights/binwidth;
  bar(centres, heights);
  set(gca, 'XLim', [-10 10]);
  
  yLim = get(gca, 'YLim')
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  axis off
  
end

if yTopStore > 1;
  yTopStore = 1;
end

for i = 1:expectedLatentDim
  
  set(ax(i, 2), 'Ylim', [0 yTopStore])
  set(ax(i, 3), 'Ylim', [0 yTopStore])
end

  
save new_demnips5.mat


