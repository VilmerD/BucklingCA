function fig = plotDesign(i_img, j_img, xx)
sizey = max(i_img);
sizex = max(j_img);
helem = sqrt(sizex*sizey/numel(xx));
fig = figure(...
    'Position', [50, 50, sizex, sizey]/helem, ...
    'MenuBar', 'none');
ax = axes(fig, ...
    'Position', [0, 0, 1, 1]);
axis(ax, 'off');
axis(ax, 'tight')
hold(ax, 'on');
dimg = sparse(i_img, j_img, xx);
imagesc(ax, dimg);
end