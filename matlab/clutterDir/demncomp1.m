clear all;
rand('seed', 1e3)
randn('seed', 1e3);
ndata = 3500;
niters = 100;
display = 2;


truenu = [2 2 2 2 2];
trueLatentDim = length(truenu);

% Tolerances
tol = 1e-5;
tautol = 1e-3;
lntautol = 1e-3;
lltol = 1e-3;

expectedLatentDim = 6;
dataDim = 6; 

trueA = [round(randn(dataDim, trueLatentDim)*40)/10];
trueBeta = 10000;

nu = 1*ones(1, expectedLatentDim);
sigma2_tau = 1*ones(1, expectedLatentDim);

for times = 1:10
   for ndata = 20:10:100
      tauSample = zeros(ndata, trueLatentDim);	
      for i = 1:trueLatentDim
         tauSample(:, i) = stgamrnd(truenu(i), 1, ndata, 1);
    end
    sSample = randn(ndata, trueLatentDim)./sqrt(tauSample);
    sNorm = nrepmat(sqrt(mean(sSample.*sSample)), 1, ndata);
    sSample = sSample./sNorm;
    
    X = sSample*trueA' + randn(ndata, dataDim)/sqrt(trueBeta);
    
    
    

    V = randn(expectedLatentDim, dataDim)*0.01;

 

    startBeta = 1/mean(diag(cov(X)));
    beta = 100*startBeta;



    nubar_tau = nrepmat(nu, 1, ndata);
    sigma2bar_tau = nrepmat(sigma2_tau, 1, ndata);

    tau = 1./sigma2bar_tau;
    lntau = zeros(ndata, expectedLatentDim);
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
	%  fprintf(1, 'Iteration %i.%i max V change: %d\n', ...
	%	  iter, counter1, maxchange)
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
	  %fprintf('Tau Iteration %i.%i, tau change %d\n', iter, taucounter, objdif);
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
      Atimes{times} = U*diag(lambda)*V';
      nutimes{times} = nu;
      betatimes{times} = beta;
      save nips1save Atimes betatimes nutimes
    end
  end
end


