%load demncompfoetal

[void, index]  = min(NU_TAU)

figure
numBins = 50;
[bins, centres] = hist(SBAR(:, index), numBins)
binwidth = centres(2)  - centres(1);

bins = bins/sum(bins);
bins = bins/binwidth;
patches = bar(centres, bins, 0.8);
xlim = get(gca, 'xlim')
x = linspace(xlim(1), xlim(2), 200);
y = tpdf(x, NU_TAU(index), SIGMA2_TAU(index));
hold on
plot(x, y);