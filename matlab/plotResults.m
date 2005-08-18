function plotResults(model, X, numOrigSignals)

% PLOTRESULTS Plot the results shown in the paper.

% GCA

if nargin < 3
    numOrigSignals = model.dataDim;
end
figNo = 1;
fontSize = 24;
fontName = 'times';
ytick = [0 0.5 1];
ylim = [0 1.2];
% NU value over which signal is considered Gaussian
nuGauss = 20;
[void, order] = sort(model.nu_tau);

% Plot the independant components, their histograms and their functional forms

figure(figNo)
figNo = figNo + 1;
nu = model.nu_tau(order);

lastIC = max(find(nu< nuGauss));
firstPC = lastIC+1;

% Orthogonalise the portion of the mixing matrix associated with the
% Gaussian directions

orthoA = model.A(:, order);
orthoA(:, end:-1:firstPC) = orthogonalise(orthoA(:, firstPC:end)')';
orthoS = (pinv(orthoA)*X')';

yTopStore = 0;
counter = 0;
orthoSLim(1) = min(min(orthoS));
orthoSLim(2) = max(max(orthoS));
plotHeight = 1/model.latentDim-0.05/model.latentDim;
figure
for i = 1:model.latentDim
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (model.latentDim - counter)/model.latentDim+0.01 0.99 plotHeight]);
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
plotHeight = 1/numOrigSignals-0.05/numOrigSignals;
figure
for i = 1:numOrigSignals
  counter = counter + 1;
  
  ax(counter, 1) = axes('position', [0.02 (numOrigSignals - counter)/numOrigSignals+0.01 0.99 plotHeight]);
  plot(1:size(X, 1), X(:, i));
  yLim = get(gca, 'YLim');
  axis off
  
    
end

centres = -4:.25:4
figure
colormap([1 1 1])
vals = hist(model.sBar(:, order(1)), centres);
vals = vals/model.numData;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
xlim = get(gca, 'xlim')
x = linspace(xlim(1), xlim(2), 200);
y = tpdf(x, nu(1), model.sigma2_tau(order(1)));
plot(x, y, 'b-');
set(gca, 'ytick', ytick)
set(gca, 'ylim', ylim)
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)

set(gca, 'ylim', [0 1.2])
figure
colormap([1 1 1])
vals = hist(model.sBar(:, order(end)), centres);
vals = vals/model.numData;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
xlim = get(gca, 'xlim')
x = linspace(xlim(1), xlim(2), 200);
y = tpdf(x, nu(end), model.sigma2_tau(order(end)));
plot(x, y, 'b-');
set(gca, 'ytick', ytick)
set(gca, 'ylim', ylim)
set(gca, 'fontname', fontName)
set(gca, 'fontsize', fontSize)
