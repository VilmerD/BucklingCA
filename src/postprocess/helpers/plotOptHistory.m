% Plots the optimization history, ie blfs

function fig = plotOptHistory(resdat)
load(resdat, 'stats', 'INPVARS')

% Make figure
fig = figure('Name', 'Obj and Const', ...
    'MenuBar', 'figure');

% Unload data
ii = 1:numel(stats);
mu = 1./[stats.lambda]';
% ii = 1:size(stats, 1);
% mu = 1./stats(:, 4:9);

% Load target blf
load(fullfile('data/compliance_reference', [INPVARS.domain, '.mat']), ...
    'lambda');
lamstar = INPVARS.sf * lambda(1);
mustar = 1/lamstar;

mu_normed = mu/mustar;

% Plot blf's and pnorm
axblfs = axes(fig);
hold(axblfs, 'on');
for k = 1:INPVARS.nevals
    plot(axblfs, ii, mu_normed(:, k), ...
        'Displayname', sprintf('\\lambda_%i', k));
end
plot(axblfs, [min(ii), max(ii)], ones(1, 2), 'k--', ...
    'DisplayName', '\lambda^*');

axis(axblfs, [min(ii), max(ii), 0, 1.2])

xlabel(axblfs, 'Iteration');
ylabel(axblfs, '$$\displaystyle\frac{\mu_j}{\mu^\star}$$', ...
    'Interpreter', 'Latex');
% lgd = legend(axblfs, ...
%     'NumColumns', 4, ...
%     'Location', 'northoutside', ...
%     'FontSize', 10);
% lgd.ItemTokenSize = 1;

fig.Position(3:4) = [400, 300];

end