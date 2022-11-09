function fig = plotDesign(i_img, j_img, xx, cmap)
if nargin < 4; cmap = 'parula'; end
sizey = max(i_img);
sizex = max(j_img);
helem = sqrt(sizex*sizey/numel(xx));
fig = figure(...
    'MenuBar', 'none', ...
    'Position', [100, 100, sizex, sizey]/helem);
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'off');
ax.Position = [0, 0, 1, 1];
colormap(ax, cmap);
axis(ax, 'tight')
dimg = sparse(i_img, j_img, xx);
imagesc(ax, dimg);
pause(1.00);
end