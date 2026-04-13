function colors = chemPlotPalette(n)
% chemPlotPalette Fixed chemistry palette for up to 16 plotted components.

palette = [ ...
    0.000 0.447 0.698
    0.835 0.369 0.000
    0.000 0.620 0.451
    0.800 0.475 0.655
    0.337 0.706 0.914
    0.769 0.306 0.322
    0.580 0.404 0.741
    0.749 0.647 0.149
    0.298 0.447 0.690
    0.867 0.518 0.322
    0.333 0.659 0.408
    0.506 0.447 0.702
    0.392 0.710 0.804
    0.855 0.545 0.765
    0.576 0.471 0.376
    0.549 0.549 0.000];

if n > size(palette, 1)
    error('chemPlotPalette supports up to 16 chemistry components.');
end

colors = palette(1:n, :);
end
