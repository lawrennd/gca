rand('seed', 1e5)
randn('seed', 1e5);
ndata = 1000;
niters = 40;
display = 2;
load \\msrc-neil01\datasets\Me'G Data\'Pierre\Spike_Set01.mat

expectedLatentDim = 25


X = EEG';
%X = X(1:100, :);

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
  V = randn(dataDim);
else
  V = randn(expectedLatentDim);
end
V = V(1:dataDim, 1:expectedLatentDim);

% Tolerances
tol = 1e-3;
tautol = 1e-3;
lntautol = 1e-3;

 
nu = 1*ones(1, expectedLatentDim);
sigma2_tau = ones(1, expectedLatentDim);

a_beta = 1e-3;
b_beta = 1e-3;
a_alpha = 1e-8;
b_alpha = 1e-3;

abar_alpha = a_alpha;
bbar_alpha = ones(1, expectedLatentDim)*b_alpha;

alpha = abar_alpha./bbar_alpha;
lnbeta = psi(abar_alpha) - log(bbar_alpha);

abar_beta = a_beta;
bbar_beta = b_beta;
beta = abar_beta./bbar_beta;
lnbeta = psi(abar_beta) - log(bbar_beta);


Sigma_V = eye(expectedLatentDim);
Sigma_V = nrepmat(Sigma_V, 3, dataDim);
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
D = eye(expectedLatentDim);
llCounter = 0;
for iter = 1:niters
  estepCount = 0;
  lntauchange = lntautol + 1;
  tauchange = tautol + 1;
  while (tauchange > tautol | lntauchange > lntautol) & estepCount < 200
    oldtau = tau;
    oldlntau = lntau;
    estepCount = estepCount + 1;
    
    maxchange = tol + 1;
    counter1 = 0;
    while maxchange > tol & counter1 < 500
      oldV = V';
      
      [sbar, Sigma_s] = bayesupdatelatent(XULinv, V, Sigma_V, tau, beta, lambda);
      [V, Sigma_V] = bayesupdateV(sbar, Sigma_s, XULinv, beta, alpha, ...
				  lambda);
      maxchange = max(max(abs(V' -oldV)));
      counter1 = counter1 + 1;
      
      if display > 1
	fprintf(1, 'V Iteration %i max change: %d\n', counter1, maxchange)
	drawnow
      end

    end
    if display > 2
      figure(3)
      for i = 1:expectedLatentDim
	for j = i+1:expectedLatentDim
	    subplot(expectedLatentDim, expectedLatentDim, (i-1)*(expectedLatentDim)+ j)
	    plot(sbar(:, i), sbar(:, j), 'rx')
	    axis image
	    axis off
	end
      end
      drawnow
    end
    [nubar_tau sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    tau = 1./sigma2bar_tau;
    for j = 1:expectedLatentDim
      lntau(:, j) = psi(nubar_tau(j)/2) ...
	  - log(nubar_tau(j)/2) ...
	  - log(sigma2bar_tau(:,j));
    end
    
  end
  
  [abar_beta bbar_beta]= bayesupdatebeta(sbar, Sigma_s, XULinv, V, Sigma_V, a_beta, ...
					 b_beta, lambda); 
  beta = abar_beta./bbar_beta;
  lnbeta = psi(abar_beta) - log(bbar_beta);
  
  
  [abar_alpha bbar_alpha] = bayesupdatealpha(V, Sigma_V, a_alpha, ...
					     b_alpha);
  
  if display
    Vstore(counter, :) = V(:)';
    betaStore(counter, :) = beta(:)';
  end
  
  % M step in latent distributions.
  [nu, sigma2_tau] = update_tauprior(nu, sigma2_tau, tau, lntau, 'fullscg');
  llCounter = llCounter + 1;
  lll(llCounter) = bicabound(V, Sigma_V, lambda, XULinv, ...
			     nu, sigma2_tau, ...
			     a_beta, b_beta, ...
			     a_alpha, b_alpha, ...
			     sbar, Sigma_s, ...
			     nubar_tau, sigma2bar_tau, ...
			     abar_beta, bbar_beta, ...
			     abar_alpha, bbar_alpha);
  if display > 1
    figure(2)
    subplot(3, 1, 1)
    plot(Vstore)
    subplot(3, 1, 2)
    plot(log10(betaStore))
    subplot(3, 1, 3)
    plot(exp(lll(1:llCounter)/ndata));
    drawnow
  end
  if llCounter > 1
    if lll(llCounter) - lll(llCounter) < 1e-6
      break
    end
  end
end





A = U*diag(lambda)*V';
figure
for i = 1:dataDim
  for j = i+1:dataDim
    subplot(dataDim, dataDim, (i-1)*(dataDim)+ j)
    plot(X(:, i), X(:, j), 'rx')
    line([ 0 A(i, i)], [0 A(j, i)])
    line([ 0 A(i, j)], [0 A(j, j)])
    axis image
  end
end
figure
for i = 1:expectedLatentDim
  for j = i+1:expectedLatentDim
    subplot(expectedLatentDim, expectedLatentDim, (i-1)*(expectedLatentDim)+ j)
    plot(sbar(:, i), sbar(:, j), 'rx')
    axis image
  end
end




figure
yTopStore = 0;
[void, index] = sort(nu);
counter = 0;
SLim(1) = min(min(sbar));
SLim(2) = max(max(sbar));
for i = index
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (expectedLatentDim - counter)/expectedLatentDim 0.45 1/expectedLatentDim]);
  %subplot(expectedLatentDim, 3, 3*counter-2);
  plot(1:ndata, sbar(:, i));
  set(gca, 'YLim', SLim)
  axis off;
  ax(counter, 2) = axes('position', [0.51 (expectedLatentDim - counter)/expectedLatentDim 0.23 1/expectedLatentDim]);

  %subplot(expectedLatentDim, 3, 3*counter-1);
  x = linspace(-4, 4, 200);
  y = tpdf(x, nu(i), sigma2_tau(i));
  plot(x, y, 'b-')
  yLim = get(gca, 'YLim');
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  set(gca, 'XTickLabel', '')%axis off
  set(gca, 'YTickLabel', '')
  %axis([-4 4 0 0.5])
  
  ax(counter, 3) = axes('position', [0.76 (expectedLatentDim - counter)/expectedLatentDim 0.23 1/expectedLatentDim]);
%  subplot(expectedLatentDim, 3, 3*counter);
  disp(nu(i))
  index = find(sbar(:, i) > -4 & sbar(:, i) < 4);
  [heights, centres] = hist(sbar(index, i), 21);
  heights = heights/sum(heights);
  binwidth = centres(2) - centres(1);
  heights = heights/binwidth;
  bar(centres, heights);
  set(gca, 'XLim', [-4 4]);
  yLim = get(gca, 'YLim');
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
  
