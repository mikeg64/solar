% Create a plot with 2 y axes using the plotyy function
figure
[ax, h1, h2] = plotyy(height/1000000, log10(temp), height/1000000, log10(dens), 'plot');

% Add title and x axis label
xlabel('Height (Mm)')
title('VALIIIc Solar Atmosphere Model')
% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'YLabel'), 'String', 'Log T')
set(get(ax(2), 'YLabel'), 'String', 'Log Density')