%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
cmap = 'parula';
source_dir = 'data/compliance_reference';
destin_dir = 'results/designs/compliance_reference';
domains = {'column', 'spire', 'twobar'};
for k = 1:numel(domains)
% Define domain and load compliance reference
domain = domains{k};
source_filename = fullfile(source_dir, sprintf('%s.mat', domain));
load(source_filename, 'sizex', 'sizey', 'xPhys');
helem = sqrt((sizex*sizey)/numel(xPhys));         % Assuming rectangular grid
switch lower(cmap)
    case {'gray', 'hot'}
        xx = 1-xPhys;
        disp('INVERTING DENSITY FIELD');
    otherwise 
        xx = xPhys;
end
% Generate geometry
genfun = str2func(sprintf('generate_%s', domain));
[X,T,i_img,j_img,solids,voids,F,freedofs] = genfun(sizex,sizey,helem,0);
% Make figure
fig = figure(...
    'Position', [50, 50, sizex, sizey]/helem, ...
    'MenuBar', 'none');
ax = axes(fig, ...
    'Position', [0, 0, 1, 1]);
hold(ax, 'on');
axis(ax, 'off');
axis(ax, 'tight');
dimg = sparse(i_img, j_img, xx);
imagesc(ax, dimg);
colormap(ax, cmap);
% Save figure
destin_filename = fullfile(destin_dir, sprintf('%s.png', domain));
saveas(fig, destin_filename);
close(fig);
end