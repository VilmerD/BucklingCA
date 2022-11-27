% Plotts the residuals

function fig = plotResiduals(jobdir)
% Load data
load(fullfile(jobdir, 'results.mat'), ...
    'res*', 'mag*');

% Plot separate figures for each 'residual'

% Stat. residual
fig = figure();
ax = axes(fig);
xlabel(ax, 'Iteration', 'Interpreter', 'LaTeX');
ylabel(ax, '$$\frac{||\textbf{Ka} - \textbf{F}||_2}{||\textbf{F}||_2}$$', ...
    'Interpreter', 'LaTeX')


relresstat = res_stat./mag_stat;
its_ca = relresstat > 1e-9;
relresstat_ca = relresstat(its_ca);
semilogy(ax, relresstat_ca, '.');

% Eigs. residual

end