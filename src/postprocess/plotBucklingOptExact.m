%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
cmap = 'parula';
farm = 'ref7';
job_nbrs = 0:11;
source_dir = fullfile('processed_data/batches/', farm);
destin_dir = 'results/designs/blfopt_exact';
for k = job_nbrs
    % Define domain and load compliance reference
    source_filename = fullfile(source_dir, sprintf('job_%i', k), 'results.mat');
    load(source_filename, 'sizex', 'sizey', 'domain', 'xPhys');
    helem = sqrt((sizex*sizey)/numel(xPhys));         % Assuming rectangular grid
    % Generate geometry
    genfun = str2func(sprintf('generate_%s', domain));
    [X,T,i_img,j_img,solids,voids,F,freedofs] = genfun(sizex,sizey,helem,0);
    % Plot design
    switch lower(cmap)
        case {'gray', 'hot'}
            xx = 1-xPhys;
            disp('INVERTING DENSITY FIELD');
        otherwise
            xx = xPhys;
    end
    fig = plotDesign(i_img, j_img, xx, cmap);
    % Save figure
    destin_filename_hf = fullfile(destin_dir, sprintf('%s%i.pdf', domain, k));
    destin_filename_lf = fullfile(destin_dir, sprintf('%s%i.png', domain, k));
    exportgraphics(fig, destin_filename_hf, 'Append', false);
    saveas(fig, destin_filename_lf);
    close(fig);
end