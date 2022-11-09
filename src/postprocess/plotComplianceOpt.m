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
load(source_filename, 'sizex', 'sizey', 'helem', 'xPhys');
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
fig = plotDesign(i_img, j_img, xx, cmap);
% Save figure
destin_filename_hf = fullfile(destin_dir, sprintf('%s.pdf', domain));
destin_filename_lf = fullfile(destin_dir, sprintf('%s.png', domain));
exportgraphics(fig, destin_filename_hf, 'Append', false);
saveas(fig, destin_filename_lf);
cla(ax);
close(fig);
clearvars -except k cmap source_dir destin_dir domains
end