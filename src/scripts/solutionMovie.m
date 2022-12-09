%% Plotting design throughout iterations
nf = size(XV, 2);
I = (T(:, 8)/(INPVARS.helem/2) + 1)/2;
J = (T(:, 7)/(INPVARS.helem/2) + 1)/2;

fig = figure();
ax = axes(fig, 'Position', [0, 0, 1, 1]);
fig.Position(3:4) = [max(J), max(I)];
mc = [0, 1];
caxis(ax, mc);

getimg = @(i) plotDesign(I, J, XV(:, i));
clear frames
for k = 1:nf
    img = getimg(k);
    cla(ax);
    imagesc(img);
    caxis(ax, mc);
    frames(k) = getframe(fig);
end

movie(frames)


function img = plotDesign(I, J, X)
% Returns an image of the design
img = reshape(X, max(I), max(J));
end

function img = plotAsymetryY(I, J, X)
% Returns an image of the asymetry in the design
Xo = reshape(X, max(I), max(J));
Xy = flip(Xo, 1);
img = (Xo - Xy).^2;
end

function img = plotAsymetryX(I, J, X)
% Returns an image of the asymetry in the design
Xo = reshape(X, max(I), max(J));
Xy = flip(Xo, 2);
img = (Xo - Xy).^2;
end