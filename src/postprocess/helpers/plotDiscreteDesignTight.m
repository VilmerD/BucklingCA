function fig = plotDiscreteDesignTight(iimg, jimg, xx, cmap)
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
dimg = sparse(iimg, jimg, xx);
imagesc(ax, dimg);
pause(0.10);

% Resize the figure
fig.Position(3:4) = [max(jimg), max(iimg)];
end