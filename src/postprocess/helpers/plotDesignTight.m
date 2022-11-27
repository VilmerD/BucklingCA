function fig = plotDesignTight(iimg, jimg, xx, cmap)
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

dimg = sparse(iimg, jimg, xx);

% Plot on axis
imagesc(ax, dimg);
pause(1.00);

% Resize the f**king figure
fig.Position(3:4) = [max(jimg), max(iimg)];
end