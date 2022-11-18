%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
cmap = 'parula';
source_dir = 'data/compliance_reference';
destin_dir = 'results/designs/compopt_ref';
domains = {'column', 'spire', 'twobar'};
for k = 1:numel(domains)
    % Define domain and load compliance reference
    domain = domains{k};
    source_filename = fullfile(source_dir, sprintf('%s.mat', domain));
    if ~isfile(source_filename); continue; end
    load(source_filename, 'INPVARS', 'xPhys');
    load(fullfile('data/domain_data/', [INPVARS.domain, '.mat']), ...
        'sizex', 'sizey');
    % Generate geometry
    genfun = str2func(sprintf('generate_%s', INPVARS.domain));
    [X,T,i_img,j_img,solids,voids,F,freedofs] = genfun(sizex,sizey,...
        INPVARS.helem,0);
    % Make figure
    fig = plotDesignTight(i_img, j_img, xPhys, cmap);
    % Save figure
    filename = fullfile(destin_dir, sprintf('%s', domain));
    exportgraphics(fig, [filename, '.pdf'], 'Append', false);
    saveas(fig, [filename, '.png']);
    close(fig);
end