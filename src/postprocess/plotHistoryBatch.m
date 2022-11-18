%% Plot the compliance and blf's
farm = 'rom15';
job_nbrs = 0:8;

source_dir = fullfile('processed_data', farm);
destin_dir = fullfile('results/history', sprintf('blfopt_%s', farm));
mkdir(destin_dir);
for k = job_nbrs
    jobdir = sprintf('job_%i', k);
    % Plot history
    source_filename = fullfile(source_dir, jobdir, 'results.mat');
    fig = plotOptHistory(source_filename);
    % Save figure
    figbase = sprintf('history_%i', k);
    filename = fullfile(destin_dir, figbase);
    exportgraphics(fig, [filename, '.pdf'], 'Append', false);
    saveas(fig, [filename, '.png']);
    close(fig);
end