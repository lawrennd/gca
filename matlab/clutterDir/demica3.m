
rand('seed', 1e4)
randn('seed', 1e4);
ndata = 100;
niters = 40;
display = 2;
truenu = [0.5 1 1.5 2];
trueLatentDim = length(truenu);
expectedLatentDim = 10;
dataDim = 10;

trueA = [round(randn(dataDim, trueLatentDim)*40)/10];
trueBeta = 1000;
for i = 1:trueLatentDim
  tauSample(:, i) = stgamrnd(truenu(i), 1, ndata, 1);
end
sSample = randn(ndata, trueLatentDim)./sqrt(tauSample);
X = sSample*trueA' + randn(ndata, dataDim)/sqrt(trueBeta);
%pre whiten data
[q, s] = eig(cov(X));
s = diag(s);
s(find(s<0)) = eps;
lambda = sqrt(s)';
sigma2_lambda = 1e-6*ones(1, dataDim);
U = q;
UT = q';
s = diag(s);
XUL = X*U*diag(1./lambda);

trueMix = (trueA'*U*diag(1./lambda))';

%X = X./nrepmat(std(X), 1, size(X, 1));

plot(XUL(:, 1), XUL(:, 2), 'rx')

% X = zero meaned data;

ndata = size(XUL, 1);
dataDim = size(XUL, 2);
%latentDim = dataDim;

%V = mixing Weights - Matrix of latentDim Rows, dataDim Columns;
if dataDim > expectedLatentDim
  V = eye(dataDim);
else
  V = eye(expectedLatentDim);
end
V = V(1:dataDim, 1:expectedLatentDim);
for j = 1:expectedLatentDim
  V(j, :) = V(j, :)/sqrt(V(j, :)*V(j, :)');
end

% Tolerances
tol = 1e-3;
tautol = 1e-3;
lntautol = 1e-3;

 
nu = 1*ones(1, expectedLatentDim);
sigma2_tau = ones(1, expectedLatentDim);

a_beta = 1e-3;
b_beta = 1e-3;

abar_beta = a_beta;
bbar_beta = b_beta*ones(1, dataDim);
beta = abar_beta./bbar_beta;
lnbeta = psi(abar_beta) - log(bbar_beta);


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
B =diag(beta);
D = eye(expectedLatentDim);
llCounter = 0;
for iter = 1:niters
  estepCount = 0;
  lntauchange = lntautol + 1;
  tauchange = tautol + 1;
  while (tauchange > tautol | lntauchange > lntautol) & estepCount < 20
    oldtau = tau;
    oldlntau = lntau;
    estepCount = estepCount + 1;
    
    maxchange = tol + 1;
    counter1 = 0;
  
    while maxchange > tol & counter1 < 100
      oldV = V';
      
      [sbar, Sigma_s] = updatelatent(XUL, V, tau, beta);
      llCounter = llCounter + 1;
      lll(llCounter) = icabound(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
				Sigma_s,  nubar_tau, sigma2bar_tau, abar_beta, bbar_beta);
      if llCounter > 1
	if lll(llCounter) < lll(llCounter-1)
	  fprintf(1, 'WARNING: bound has decreased by %d after S update\n', lll(llCounter-1) ...
		    - lll(llCounter))
	end
      end
      V = gsupdateV(sbar, Sigma_s, XUL);
      llCounter = llCounter + 1;
      lll(llCounter) = icabound(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
				Sigma_s,  nubar_tau, sigma2bar_tau, abar_beta, bbar_beta);
      if lll(llCounter) < lll(llCounter-1)
	fprintf(1, 'WARNING: bound has decreased by %d after V update\n', lll(llCounter-1) ...
		- lll(llCounter))
      end
      
      maxchange = max(max(abs(V' -oldV)));
      counter1 = counter1 + 1;
      
      if display > 1
	fprintf(1, 'V Iteration %i max change: %d\n', counter1, ...
		maxchange)
	fprintf(1, 'Likelihood bound: %d\n', lll(llCounter));
	drawnow
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
    end
    
    [abar_beta bbar_beta]= updatebeta(sbar, Sigma_s, XUL, V, a_beta, ...
				      b_beta, lambda); 
    beta = abar_beta./bbar_beta;
    lnbeta = psi(abar_beta) - log(bbar_beta);
    
    llCounter = llCounter + 1;
    lll(llCounter) = icabound(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
			      Sigma_s,  nubar_tau, sigma2bar_tau, abar_beta, bbar_beta);
    if lll(llCounter) < lll(llCounter-1)
      fprintf(1, 'WARNING: bound has decreased by %d after BETA update\n', lll(llCounter-1) ...
	      - lll(llCounter));
    end
    
    [nubar_tau sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    tau = 1./sigma2bar_tau;
    for j = 1:expectedLatentDim
      lntau(:, j) = psi(nubar_tau(j)/2) ...
	  - log(nubar_tau(j)/2) ...
	  - log(sigma2bar_tau(:,j));
    end

    llCounter = llCounter + 1;
    lll(llCounter) = icabound(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
			      Sigma_s,  nubar_tau, sigma2bar_tau, abar_beta, bbar_beta);
    if lll(llCounter) < lll(llCounter-1)
      fprintf(1, 'WARNING: bound has decreased by %d after TAU update\n', lll(llCounter-1) ...
	      - lll(llCounter))
    end
    tauchange = max(max(abs(tau - oldtau)));
    lntauchange = max(max(abs(lntau - oldlntau)));
    fprintf(1, 'tau Iteration %i tau change: %d lntau change %d\n', ...
	    estepCount, tauchange, lntauchange)
    fprintf(1, 'Likelihood bound: %d\n', lll(llCounter));
    
    counter = counter + 1;
    if display
      Vstore(counter, :) = V(:)';
      betaStore(counter, :) = beta(:)';
    end
  end
  
  % M step in latent distributions.
  [nu, sigma2_tau] = update_tauprior(nu, sigma2_tau, tau, lntau, 'fullscg');
  llCounter = llCounter + 1;
  lll(llCounter) = icabound(V, XUL, nu, sigma2_tau, a_beta, b_beta, sbar, ...
			    Sigma_s,  nubar_tau, sigma2bar_tau, abar_beta, bbar_beta);
  if lll(llCounter) < lll(llCounter-1)
    fprintf(1, 'WARNING: bound has decreased by %d after latent parameter update\n', lll(llCounter-1) ...
	    - lll(llCounter))
  end
  
  if display > 1
    figure(2)
    subplot(3, 1, 1)
    plot(Vstore)
    subplot(3, 1, 2)
    plot(log10(betaStore))
    subplot(3, 1, 3)
    plot(exp(lll(llCounter)/ndata));
    drawnow
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
  
