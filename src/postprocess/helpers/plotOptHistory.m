% Plots the optimization history, ie blfs

function fig = plotOptHistory(resdat)
load(resdat, 'stats', 'INPVARS')

% Make figure
fig = figure('Name', 'Obj and Const', ...
    'MenuBar', 'figure');

% Unload data
ii = 1:numel(stats);
blfs = 1./[stats.lambda]';

% Load target blf
load(fullfile('data/compliance_reference', [INPVARS.domain, '.mat']), ...
    'lambda');
lamstar = INPVARS.sf * lambda(1);
mustar = 1/lamstar;

% Plot blf's and pnorm
axblfs = axes(fig);
hold(axblfs, 'on');
for k = 1:INPVARS.nevals
    plot(axblfs, ii, blfs(:, k), ...
        'Displayname', sprintf('\\lambda_%i', k));
end
plot(axblfs, [min(ii), max(ii)], mustar*ones(1, 2), 'k--', ...
    'DisplayName', '\lambda^*');

axis(axblfs, [min(ii), max(ii), 0, mustar*1.2])

xlabel(axblfs, 'Iteration');
ylabel(axblfs, 'BLFs');
% lgd = legend(axblfs, ...
%     'NumColumns', 4, ...
%     'Location', 'northoutside', ...
%     'FontSize', 10);
% lgd.ItemTokenSize = 1;

fig.Position(3:4) = [400, 300];

end