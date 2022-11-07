%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
cmap = 'parula';
source_dir = 'processed_data/ref2/';
job_nbrs = 6:8;
destin_dir = 'results/designs/blfopt_exact';
domains = {'column', 'spire', 'twobar'};
for k = 1:numel(domains)
    % Define domain and load compliance reference
    domain = domains{k};
    jobname = sprintf('job_%i', job_nbrs(k));
    source_filename = fullfile(source_dir, jobname, 'results.mat');
    load(source_filename, 'xPhys');
    switch lower(domain)
        case 'column'
            sizex = 240;
            sizey = 120;
        case 'spire'
            sizex = 240;
            sizey = 120;
        case 'twobar'
            sizex = 120;
            sizey = 280;
    end
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
    fig = plotDesign(i_img, j_img, xx);
    figure(fig);
    ax = gca();
    colormap(ax, cmap);
    % Save figure
    destin_filename = fullfile(destin_dir, sprintf('%s.png', domain));
    saveas(fig, destin_filename);
    close(fig);
end