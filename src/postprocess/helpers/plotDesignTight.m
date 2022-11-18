function fig = plotDesignTight(i_img, j_img, xx, cmap)
% PLOTDESIGNTIGHT(i_img, j_img, xx) plots the design specified in xx with
% coordinates i_img and j_img

% Default cmap and
if nargin < 4; cmap = 'parula'; end

% Hot or gray cmap requires inverting the field
if sum(strcmpi(cmap, {'hot', 'gray'})) > 0 
    disp('INVERTING DENSITY FIELD');
    xx = 1-xx;
end
% Make figure of correct size
sizey = max(i_img);
sizex = max(j_img);
helem = sqrt(sizex*sizey/numel(xx));
fig = figure(...
    'MenuBar', 'none', ...
    'Position', [100, 100, sizex, sizey]/helem);

% Make axis in figure
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'off');
ax.Position = [0, 0, 1, 1];
colormap(ax, cmap);
axis(ax, 'tight')

% Make image
dimg = sparse(i_img, j_img, xx);

% Plot on axis
imagesc(ax, dimg);
pause(0.50);
end