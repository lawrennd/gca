load demncompfoetal

[void, index]  = min(nu)

figure
numBins = 50;
[bins, centres] = hist(sbar(:, index), numBins)
binwidth = centres(2)  - centres(1);

bins = bins/sum(bins);
bins = bins/binwidth;
patches = bar(centres, bins, 0.8);
xlim = get(gca, 'xlim')
x = linspace(xlim(1), xlim(2), 200);
y = tpdf(x, nu(index), sigma2_tau(index));
hold on
plot(x, y);