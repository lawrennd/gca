rand('seed', 1e5)
randn('seed', 1e5);
ndata = 1000;
niters = 40;
display = 2;
load c:\datasets\Me'G Data\'Pierre\Spike_Set01.mat

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
XUL = X*U*diag(1./lambda);


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
      V = updateV(sbar, Sigma_s, XUL);
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
for i = index
  counter = counter + 1;
  subplot(expectedLatentDim, 3, 3*counter-2);
  plot(1:ndata, sbar(:, i));
  axis off;
  subplot(expectedLatentDim, 3, 3*counter-1);
  x = linspace(-4, 4, 200);
  y = tpdf(x, nu(i), sigma2_tau(i));
  plot(x, y, 'b-')
  yLim = get(gca, 'YLim');
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  
  %axis([-4 4 0 0.5])
  
  subplot(expectedLatentDim, 3, 3*counter);
  disp(nu(i))
  [heights, centres] = hist(sbar(:, i), 20);
  heights = heights/sum(heights);
  bar(centres, heights);
  %axis([-4 4 0 0.5])
  yLim = get(gca, 'YLim');
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  
end

for i = 1:expectedLatentDim
  subplot(expectedLatentDim, 3, 3*i-1);
  set(gca, 'Ylim', [0 yTopStore])
  subplot(expectedLatentDim, 3, 3*i);
  set(gca, 'Ylim', [0 yTopStore])
end
  
