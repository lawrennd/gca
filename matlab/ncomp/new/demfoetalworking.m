clear all;
HOME = getenv('HOME');
addpath([HOME '/mlprojects/gca/matlab'])
rand('seed', 1e3)
randn('seed', 1e3);
display = 2;
load 'c:\datasets\ICA\ECG Data\foetal_ecg.dat'

global NDATA
global DATADIM
global LATENTDIM

global X

global BETA
global NU_TAU
global SIGMA2_TAU
global LAGRANGE % Lagrange multipliers

global NUBAR_TAU
global SIGMA2BAR_TAU

global TAU
global LNTAU

global A
global SBAR
global SIGMA_S

LATENTDIM = 12;
X = foetal_ecg(:, 2:9);%(1:10:end, :);%(:, 2:9); %6500)';
NDATA = size(X, 1);
DATADIM = size(X, 2);
% X = zero meaned data;
meanX = mean(X, 1);
for i = 1:NDATA
  X(i, :)  = X(i, :) - meanX;
end
for i = 1:DATADIM
  X(:, i) = X(:, i)/std(X(:, i));
end

%pre whiten data
%[q, s] = eig(cov(X));
%s = diag(s);
%s(find(s<0)) = eps;
%lambda = sqrt(s)';
%U = q;
%UT = q';
%s = diag(s);
A = randn( DATADIM, LATENTDIM)*0.1;

plot(X(:, 1), X(:, 2), 'rx')

% X = zero meaned data;

%NDATA = size(X, 1);
%DATADIM = size(X, 2);
%LATENTDIM = DATADIM;


% Tolerances
tol = 1e-5;
tautol = 1e-3;
Lntautol = 1e-3;
lltol = 1e-3;
 
NU_TAU = 3*ones(1, LATENTDIM);
SIGMA2_TAU = 1/3*ones(1, LATENTDIM);

BETA = min(eig(cov(X)));



NUBAR_TAU = nrepmat(NU_TAU, 1, NDATA);
SIGMA2BAR_TAU = nrepmat(SIGMA2_TAU, 1, NDATA);

TAU = 1./SIGMA2BAR_TAU;
for j = 1:LATENTDIM
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
      - log(NUBAR_TAU(j)/2) ...
      - log(SIGMA2BAR_TAU(:,j));
end
min_tau = 2.001;
counter = 1;
Astore = A(:)';
betaStore = BETA(:)';
nu_tauStore = NU_TAU(:)';
D = eye(LATENTDIM);
counter = 1;
SBAR = zeros(NDATA, LATENTDIM);
updatelatent
stupdatetau
TAU = 1./SIGMA2BAR_TAU;
for j = 1:LATENTDIM
  LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
      - log(NUBAR_TAU(j)/2) ...
      - log(SIGMA2BAR_TAU(:,j));
end
LAGRANGE = zeros(size(NU_TAU));
lll(1) = sticabound/NDATA;
llldiff = 1;
order = randperm(5);
while  counter < 400
  for i = order
    switch i
     case 1
      updatelatent
     case 2
      stupdatetau
      TAU = 1./SIGMA2BAR_TAU;
      for j = 1:LATENTDIM
	LNTAU(:, j) = psi(NUBAR_TAU(j)/2) ...
	    - log(NUBAR_TAU(j)/2) ...
	    - log(SIGMA2BAR_TAU(:,j));
      end
     case 3
      updateA
     case 4
      updatebeta
     case 5 
      stupdatetauprior('scg', min_tau)
    end
  end
    counter = counter + 1;
    lll(counter) = sticabound/NDATA;
    llldiff = lll(counter) - lll(counter-1);
    fprintf('Iteration %i.%i, log likelihood change: %d\n', counter, i,llldiff)
    
  if display
    Astore(counter, :) = A(:)';
    betaStore(counter, :) = BETA(:)';
    nu_tauStore(counter, :) = NU_TAU(:)';
  end
  
  if display > 1
    figure(1)
    subplot(4, 1, 1)
    plot(Astore)
    subplot(4, 1, 2)
    plot(log10(betaStore))
    subplot(4, 1, 3)
    plot(log10(nu_tauStore))
    subplot(4, 1, 4)
    plot(1:counter, lll(1:counter));
    drawnow
  end  
  save ncompfoetalsavequick A BETA NU_TAU
%  if abs(llldiff) < 1e-3
%    break
%  end
  pause(0.1)
end

figure(3)
for i = 1:LATENTDIM
  for j = i+1:LATENTDIM
    subplot(LATENTDIM, LATENTDIM, (i-1)*(LATENTDIM)+ j)
    plot(SBAR(:, i), SBAR(:, j), 'rx')
    axis image
  end
end


figure(4)
clf
yTopStore = 0;
[void, index] = sort(NU_TAU);
counter = 0;
SBARLim(1) = min(min(SBAR));
SBARLim(2) = max(max(SBAR));
plotHeight = 1/LATENTDIM-0.05/LATENTDIM;
for i = index
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (LATENTDIM - counter)/LATENTDIM+0.01 0.77 plotHeight]);
  plot(1:size(SBAR, 1), SBAR(:, i));
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
  ax(counter, 2) = axes('position', [0.81 (LATENTDIM - counter)/LATENTDIM+0.01 0.08 plotHeight]);
  x = linspace(-4, 4, 200);
%  y = tpdf(x, NU_TAU(i), SIGMA2_TAU(i));
%  plot(x, y, 'b-')
  yLim = get(gca, 'YLim');
  if yLim(2) > yTopStore
    yTopStore = yLim(2);
  end
  set(gca, 'XLim', [-10 10]);
  set(gca, 'XTickLabel', '')%axis off
  set(gca, 'YTickLabel', '')
  ax(counter, 3) = axes('position', [0.91 (LATENTDIM - counter)/LATENTDIM+0.01 0.08 plotHeight]);
  disp(NU_TAU(i))
  index2 = find(SBAR(:, i) > -10 & SBAR(:, i) < 10);
  index2 = 1:length(index2);
  [heights, centres] = hist(SBAR(index2, i), 30);
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

for i = 1:LATENTDIM
  
  set(ax(i, 2), 'Ylim', [0 yTopStore])
  set(ax(i, 3), 'Ylim', [0 yTopStore])
end

  
save demncompfoetalquick.mat


