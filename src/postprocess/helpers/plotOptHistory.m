% Plots the optimization history, ie blfs

function fig = plotOptHistory(resdat)
load(resdat, 'stats', 'INPVARS')

% Make figure
fig = figure('Name', 'Obj and Const', ...
    'MenuBar', 'figure');

% Unload data
ii = 1:size(stats, 1);
blfs = stats(ii, 4+(0:INPVARS.nevals-1));

% Load target blf
load(fullfile('data/compliance_reference', [INPVARS.domain, '.mat']), ...
    'lambda');
lamstar = INPVARS.sf * lambda(1);

% Plot blf's and pnorm
axblfs = axes(fig);
hold(axblfs, 'on');
for k = 1:INPVARS.nevals
    plot(axblfs, ii, blfs(:, k), ...
        'Displayname', sprintf('\\lambda_%i', k));
end
plot(axblfs, [min(ii), max(ii)], lamstar*ones(1, 2), 'k--', ...
    'DisplayName', '\lambda^*');

axis(axblfs, [min(ii), max(ii), 0, mean(blfs(:, end))*1.2])

xlabel(axblfs, 'Iteration');
ylabel(axblfs, 'BLFs');
% lgd = legend(axblfs, ...
%     'NumColumns', 4, ...
%     'Location', 'northoutside', ...
%     'FontSize', 10);
% lgd.ItemTokenSize = 1;

fig.Position(3:4) = [300, 300];

end

function fig = plotOptHistoryAlt(resdat)
load(resdat, 'stats', 'INPVARS')

% Make figure
fig = figure('Name', 'Obj and Const', ...
    'MenuBar', 'figure');
t = tiledlayout(fig, 2, 3, ...
    'TileSpacing', 'compact');

% Unload data
ii = 1:size(stats, 1);
g0 = stats(ii, 1);
g11 = stats(ii, 2);
g12 = stats(ii, 3);
blfs = stats(ii, 4+(0:INPVARS.nevals-1));

% Load target blf
load(fullfile('data/compliance_reference', [INPVARS.domain, '.mat']), ...
    'lambda');
lamstar = INPVARS.sf * lambda(1);

% Plot objective
axobj = nexttile(t, 1, [1, 1]);
plot(axobj, ii, g0);

xlabel(axobj, 'Iteration');
ylabel(axobj, 'Objective');

% Plot blf's and pnorm
axblfs = nexttile(t, 2, [2, 2]);
hold(axblfs, 'on');
for k = 1:INPVARS.nevals
    plot(axblfs, ii, blfs(:, k), ...
        'Displayname', sprintf('\\lambda_%i', k));
end
plot(axblfs, [min(ii), max(ii)], lamstar*ones(1, 2), 'k--', ...
    'DisplayName', '\lambda^*');

axis(axblfs, [min(ii), max(ii), 0, mean(blfs(:, end))*1.2])

xlabel(axblfs, 'Iteration');
ylabel(axblfs, 'BLFs');
legend(axblfs, ...
    'NumColumns', 2, ...
    'Location', 'southwest');

% Plot Constraints
axcons = nexttile(t, 4, [1, 1]);
hold(axcons, 'on');
plot(axcons, g12, 'Displayname', 'BLF');
plot(axcons, g11, 'Displayname', 'V');

axis(axcons, [min(ii), max(ii), -0.20, 0.20])

xlabel(axcons, 'Iteration');
ylabel(axcons, 'Constraint');
legend(axcons);
end