fontSize = 24;
fontName = 'times';
ytick = [0 0.5 1];
ylim = [0 1.2];
% NU value over which signal is considered Gaussian
nuGauss = 20;
[void, order] = sort(NU_TAU);

% Plot the independant componets, their histograms and their functional forms

figure(figNo)
figNo = figNo + 1;
nu = NU_TAU(order);

lastIC = max(find(nu< nuGauss));
firstPC = lastIC+1;

% Orthogonalise the portion of the mixing matrix associated with the
% Gaussian directions
orthoA = A(:, order);
orthoA(:, end:-1:firstPC) = orthogonalise(orthoA(:, firstPC:end)')';
orthoS = (pinv(orthoA)*X')';

yTopStore = 0;
counter = 0;
orthoSLim(1) = min(min(orthoS));
orthoSLim(2) = max(max(orthoS));
plotHeight = 1/LATENTDIM-0.05/LATENTDIM;
figure
for i = 1:LATENTDIM
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (LATENTDIM - counter)/LATENTDIM+0.01 0.99 plotHeight]);
  plot(1:size(orthoS, 1), orthoS(:, i));
  yLim = get(gca, 'YLim');
  
  set(gca, 'YLim', [-10 10]);
  set(gca, 'XTickLabel', '')
  set(gca, 'YTick', [yLim(1) 0 yLim(2)])
  axis off
end

yTopStore = 0;
counter = 0;
XLim(1) = min(min(X));
XLim(2) = max(max(X));
plotHeight = 1/DATADIM-0.05/DATADIM;
figure
for i = 1:DATADIM
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (DATADIM - counter)/DATADIM+0.01 0.99 plotHeight]);
  plot(1:size(X, 1), X(:, i));
  yLim = get(gca, 'YLim');
  axis off
  
%  set(gca, 'YLim', [-100 100]);
%  set(gca, 'XTickLabel', '')
%  set(gca, 'YTick', [yLim(1) 0 yLim(2)])
%  axis 
    
end

centres = -4:.25:4
figure
colormap([1 1 1])
vals = hist(SBAR(:, order(1)), centres);
vals = vals/NDATA;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
y = tpdf(x, nu(1), SIGMA2_TAU(order(1)));
plot(x, y, 'b-');
set(gca, 'ytick', ytick)
set(gca, 'ylim', ylim)
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)

set(gca, 'ylim', [0 1.2])
figure
colormap([1 1 1])
vals = hist(SBAR(:, order(end)), centres);
vals = vals/NDATA;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
y = tpdf(x, nu(end), SIGMA2_TAU(order(end)));
plot(x, y, 'b-');
set(gca, 'ytick', ytick)
set(gca, 'ylim', ylim)
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
