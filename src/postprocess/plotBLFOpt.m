%% Plot the compliance optimized shapes
% Assuming rectangular grid
% Define some parameters
farm = 'ref022';

cmap = 'parula';
plotstyle = 'continuous';
source_dir = fullfile('processed_data/', farm);
destin_dir = fullfile('results/designs', sprintf('blfopt_%s', farm));
mkdir(destin_dir);

nj = numel(dir(fullfile(source_dir, 'job_*')));
for k = 0:nj-1
    % Define domain and load compliance reference
    source_filename = fullfile(source_dir, sprintf('job_%i', k), 'results.mat');
    if ~isfile(source_filename); continue; end
    load(source_filename, 'INPVARS', 'x', 'xPhys');
    load(fullfile('data/domain_data/', [INPVARS.domain, '.mat']), ...
        'sizex', 'sizey');
    % Generate element connectivity and stuff for plot
    genfun = str2func(sprintf('generate_%s', INPVARS.domain));
    [X,T,i_img,j_img,solids,voids,F,freedofs] = genfun(sizex,sizey,...
        INPVARS.helem,0);
    % Plot design
    if strcmpi(plotstyle, 'continuous')
        [F, edofMat] = makeFilter(INPVARS.rmin, INPVARS.xi, T, INPVARS.helem, solids, voids);
        xtilde = F(x);
        xphys = heaviside_threshold(xtilde, INPVARS.betaend, 0.50);
        edens = xphys(edofMat);
        fig = plotContinuousDesignTight(i_img, j_img, edens, cmap);
    elseif strcmpi(plotstyle, 'discrete')
        fig = plotDiscreteDesignTight(i_img, j_img, xPhys, cmap);
    end
    % Save figure
    filename = fullfile(destin_dir, sprintf('%s%i', INPVARS.domain, k));
    exportgraphics(fig, [filename, '.pdf'], 'Append', false);
    saveas(fig, [filename, '.png']);
    close(fig);
end