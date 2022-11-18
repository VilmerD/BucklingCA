%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
farm = 'rom15';
job_nbrs = 0:8;

cmap = 'parula';
source_dir = fullfile('processed_data/', farm);
destin_dir = fullfile('results/designs', sprintf('blfopt_%s', farm));
mkdir(destin_dir);
for k = job_nbrs
    % Define domain and load compliance reference
    source_filename = fullfile(source_dir, sprintf('job_%i', k), 'results.mat');
    if ~isfile(source_filename); continue; end
    load(source_filename, 'INPVARS', 'xPhys');
    load(fullfile('data/domain_data/', [INPVARS.domain, '.mat']), ...
        'sizex', 'sizey');
    % Generate geometry
    genfun = str2func(sprintf('generate_%s', INPVARS.domain));
    [X,T,i_img,j_img,solids,voids,F,freedofs] = genfun(sizex,sizey,...
        INPVARS.helem,0);
    % Plot design
    fig = plotDesignTight(i_img, j_img, xPhys, cmap);
    % Save figure
    filename = fullfile(destin_dir, sprintf('%s%i', INPVARS.domain, k));
    exportgraphics(fig, [filename, '.pdf'], 'Append', false);
    saveas(fig, [filename, '.png']);
    close(fig);
end