% NU value over which signal is considered Gaussian
nuGauss = 20;
[void, order] = sort(NU_TAU);

figNo = 2;
if LATENTDIM < 11
  figure(figNo)
  figNo = figNo + 1;
  % Plot of all the diferent Indepenent components vs one another.
  for i = 1:LATENTDIM
    for j = i+1:LATENTDIM
      subplot(LATENTDIM, LATENTDIM, (i-1)*(LATENTDIM)+ j)
      plot(SBAR(:, order(i)), SBAR(:, order(j)), 'rx')
      axis image
    end
  end
end

% Plot the independant componets, their histograms and their functional forms
figure(figNo)
figNo = figNo + 1;
clf
yTopStore = 0;
[void, order] = sort(NU_TAU);
counter = 0;
SBARLim(1) = min(min(SBAR));
SBARLim(2) = max(max(SBAR));
plotHeight = 1/LATENTDIM-0.05/LATENTDIM;
for i = order
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

  


if any(NU_TAU>nuGauss)
  figure(figNo)
  figNo = figNo + 1;
  nu = NU_TAU(order);
  
  lastIC = max(find(nu< nuGauss));
  firstPC = lastIC+1;
  
  % Orthogonalise the portion of the mixing matrix associated with the
  % Gaussian directions
  PCS = find(NU_TAU> nuGauss);
  orthoA = A;
  orthoA(:, PCS) = orthogonalise(orthoA(:, PCS)')';
  orthoS = (pinv(orthoA)*X')';
  
  yTopStore = 0;
  counter = 0;
  orthoSLim(1) = min(min(orthoS));
  orthoSLim(2) = max(max(orthoS));
  plotHeight = 1/LATENTDIM-0.05/LATENTDIM;
  for i = 1:LATENTDIM
    counter = counter + 1;
    
    ax(counter, 1) = axes('position', [0.02 (LATENTDIM - counter)/LATENTDIM+0.01 0.77 plotHeight]);
    plot(1:size(orthoS, 1), orthoS(:, order(i)));
    yLim = get(gca, 'YLim');
    
    set(gca, 'YLim', [-10 10]);
    set(gca, 'XTickLabel', '')
    set(gca, 'YTick', [yLim(1) 0 yLim(2)])
    axis 
    ax(counter, 2) = axes('position', [0.81 (LATENTDIM - counter)/LATENTDIM+0.01 0.08 plotHeight]);
    x = linspace(-4, 4, 200);
    %  y = tpdf(x, nu(i), SIGMA2_TAU(order(i)));
    %  plot(x, y, 'b-')
    yLim = get(gca, 'YLim');
    if yLim(2) > yTopStore
      yTopStore = yLim(2);
    end
    set(gca, 'XLim', [-10 10]);
    set(gca, 'XTickLabel', '')%axis off
    set(gca, 'YTickLabel', '')
    ax(counter, 3) = axes('position', [0.91 (LATENTDIM - counter)/LATENTDIM+0.01 0.08 plotHeight]);
    disp(nu(i))
    index2 = find(orthoS(:, order(i)) > -10 & orthoS(:, order(i)) < 10);
    index2 = 1:length(index2);
    [heights, centres] = hist(orthoS(index2, order(i)), 30);
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
end
centres = -4:.25:4
figure
colormap([1 1 1])
vals = hist(SBAR(:, order(1)), centres);
vals = vals/NDATA;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
y = tpdf(x, NU_TAU(order(1)), SIGMA2_TAU(order(1)));
plot(x, y, 'b-');

set(gca, 'ylim', [0 1.2])
figure
colormap([1 1 1])
vals = hist(SBAR(:, order(end)), centres);
vals = vals/NDATA;
vals = vals/(centres(2) - centres(1))
bar(centres(2:end-1), vals(2:end-1));
hold on
y = tpdf(x, NU_TAU(order(end)), SIGMA2_TAU(order(end)));
plot(x, y, 'b-');
%set(gca, 'ytick', [0 0.1 0.2 0.3 0.4])
set(gca, 'ylim', [0 1.2])
