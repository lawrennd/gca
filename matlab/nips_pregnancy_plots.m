load nips_pregnancy

dataDim = size(X, 2);
figure(1)
for i = 1:dataDim
  ax(i) = subplot(dataDim, 1, i)
  plot(1:ndata, X(:, i));
  yLim = get(ax(i), 'YLim')
  if yLim(1) > -5
    yLim(1) = -5;
  end
  if yLim(2) < 5
    yLim(2) = 5;
  end
  set(ax(i), 'YLim', yLim);
  set(ax(i), 'YTick', [yLim(1) 0 yLim(2)])
  yTickLabel = get(ax(i), 'YTickLabel');
  yTickLabel(1, :) = char(zeros(1, length(yTickLabel(1, :))));
  set(ax(i), 'YTickLabel', yTickLabel)
  set(ax(i), 'XTick', [])
end
printlatex('z:\tex\projects\ica\diagrams\pregnancy_data.eepic', 12, 'top', 'y')  

% First print raw solutions
[nu, order] = sort(nu);

A = A(:, order);

s = (pinv(A)*X')';
ndata = size(s, 1);
latentDim = size(s, 2);


figure(2)
for j = 1:latentDim
  ax(j) = subplot(latentDim, 1, j)
  plot(1:ndata, s(:, j));
  yLim = get(ax(j), 'YLim')
  if yLim(1) > -5
    yLim(1) = -5;
  end
  if yLim(2) < 5
    yLim(2) = 5;
  end
  set(ax(j), 'YLim', yLim);
  set(ax(j), 'YTick', [yLim(1) 0 yLim(2)])
  yTickLabel = get(ax(j), 'YTickLabel');
  yTickLabel(1, :) = char(zeros(1, length(yTickLabel(1, :))));
  set(ax(j), 'YTickLabel', yTickLabel)
  set(ax(j), 'XTick', [])
end
printlatex('z:\tex\projects\ica\diagrams\pregnancy_raw.eepic', 12, 'top', 'y')  

lastIC = max(find(nu< 100));
firstPC = lastIC+1;
orthoA = A;
orthoA(:, end:-1:firstPC) = orthogonalise(A(:, firstPC:end)')';

orthoS = (pinv(orthoA)*X')';
ndata = size(s, 1);
latentDim = size(s, 2);


figure(3)
for j = 1:latentDim
  ax(j) = subplot(latentDim, 1, j)
  plot(1:ndata, orthoS(:, j));
  yLim = get(ax(j), 'YLim')
  if yLim(1) > -5
    yLim(1) = -5;
  end
  if yLim(2) < 5
    yLim(2) = 5;
  end
  set(ax(j), 'YLim', yLim);
  set(ax(j), 'YTick', [yLim(1) 0 yLim(2)])
  yTickLabel = get(ax(j), 'YTickLabel');
  yTickLabel(1, :) = char(zeros(1, length(yTickLabel(1, :))));
  set(ax(j), 'YTickLabel', yTickLabel)
  set(ax(j), 'XTick', [])
end
  

printlatex('z:\tex\projects\ica\diagrams\pregnancy_subspace.eepic', 12, 'top', 'y')  

figure(4)
x = linspace(-5, 5, 50);
y = tpdf(x, nu(1), 1);
plot(x, y, 'b-')
set(gca, 'xtick', [-5 0 5])
set(gca, 'ytick', [0 0.2 0.4])
set(gca, 'ylim', [0 0.4])
hold on
%printlatex('z:\tex\projects\ica\diagrams\pregnancy_IC1_distribution.eepic', 5, 'top', 'y')  

figure(4)
x = linspace(-5, 5, 50);
y = tpdf(x, nu(6), 1);
plot(x, y, 'b:')
set(gca, 'xtick', [-5 0 5])
set(gca, 'ytick', [0 0.2 0.4])
printlatex('z:\tex\projects\ica\diagrams\pregnancy_IC1PC_distribution.eepic', 3, 'top', 'y')  

