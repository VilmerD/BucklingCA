function fig = plotContinuousDesignTight(iimg, jimg, xx, cmap)
% PLOTDESIGNTIGHT(i_img, j_img, xx) plots the design specified in xx with
% coordinates i_img and j_img

% Default cmap and
if nargin < 4; cmap = 'parula'; end

% Hot or gray cmap requires inverting the field
if sum(strcmpi(cmap, {'hot', 'gray'})) > 0 
    disp('INVERTING DENSITY FIELD');
    xx = 1-xx;
end
% Make figure
fig = figure(...
    'MenuBar', 'none');

% Make axis in figure
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'off');
ax.Position = [0, 0, 1, 1];
colormap(ax, cmap);
axis(ax, 'tight')

% Prepare figure
% Prepare element coordiantes
ex = jimg + [-0.50, -0.50,  0.50,  0.50];
ey = iimg + [-0.50,  0.50,  0.50, -0.50];
% Sort density values such that the biggest one is first, makes the plot
% better looking
for k = 1:size(xx, 1)
    [~, I] = min(xx(k, :));
    xx(k, :) = circshift(xx(k, :), 1-I);
    ex(k, :) = circshift(ex(k, :), 1-I);
    ey(k, :) = circshift(ey(k, :), 1-I);
end
patch(ax, ex', ey', xx', 'LineStyle', 'none');
pause(0.10);

% Resize the figure
fig.Position(3:4) = [max(jimg), max(iimg)];
end