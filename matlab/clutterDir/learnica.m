function [A, beta, nu, sbar] = learnica(X, latentDim);

ndata = size(X, 1);
dataDim = size(X, 2);
niters = 100;
display = 2;

% Tolerances
tol = 1e-5;
tautol = 1e-3;
lntautol = 1e-3;
lltol = 1e-3;

% X zero mean data;
meanX = mean(X, 1);
for i = 1:ndata
  X(i, :)  = X(i, :) - meanX;
end

% pre whiten data
[q, s] = eig(cov(X));
s = diag(s);
s(find(s<0)) = eps;
lambda = sqrt(s)';
sigma2_lambda = 1e-6*ones(1, dataDim);
U = q;
UT = q';
s = diag(s);
XULinv = X*U*diag(1./lambda);

V = randn(latentDim, dataDim)*0.01; 

startBeta = 1/mean(diag(cov(X)));
beta = 100*startBeta;

nu = 1*ones(1, latentDim);
sigma2_tau = 1*ones(1, latentDim);

nubar_tau = nrepmat(nu, 1, ndata);
sigma2bar_tau = nrepmat(sigma2_tau, 1, ndata);

tau = 1./sigma2bar_tau;
lntau = zeros(ndata, latentDim);
for j = 1:latentDim
  lntau(:, j) = psi(nubar_tau(j)/2) ...
      - log(nubar_tau(j)/2) ...
      - log(sigma2bar_tau(:,j));
end

counter = 1;
Vstore = V(:)';
betaStore = beta(:)';
nuStore = nu(:)';
D = eye(latentDim);
llCounter = 1;
sbar = zeros(ndata, latentDim);
[sbar, Sigma_s] = updatelatent(XULinv, V, tau, beta, lambda);
lll(llCounter) = sticabound(V, lambda, XULinv, beta, ...
			    nu, sigma2_tau, ...
			    sbar, Sigma_s, ...
			    nubar_tau, sigma2bar_tau);
for iter = 1:niters
  estepCount = 0;
  lntauchange = lntautol + 1;
  tauchange = tautol + 1;
  
  maxchange = tol + 1;
  counter1 = 0;
  
  while maxchange > tol & counter1 < 500
    
    oldV = V';
    [sbar, Sigma_s] = updatelatent(XULinv, V, tau, beta, ...
				   lambda);
    V = updateV(sbar, Sigma_s, XULinv);
    [nubar_tau sigma2bar_tau] = stupdatetau(sbar, Sigma_s, nu, sigma2_tau);
    tau = 1./sigma2bar_tau;
    for j = 1:latentDim
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
    %    % End of SV loop
  end
  
  beta = updatebeta(sbar, Sigma_s, XULinv, V, lambda); 
  
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
  for j = 1:latentDim
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
    for j = 1:latentDim
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
  %if llldiff < tol
  %  break
  %end
end
A =  U*diag(lambda)*V';

