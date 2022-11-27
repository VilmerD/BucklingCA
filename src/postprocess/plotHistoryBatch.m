%% Plot the compliance and blf's
farm = 'ref010';

source_dir = fullfile('processed_data', farm);
destin_dir = fullfile('results/history', sprintf('blfopt_%s', farm));
mkdir(destin_dir);

nj = numel(dir(fullfile(source_dir, 'job_*')));
for k = 0:nj-1
    jobdir = sprintf('job_%i', k);
    % Plot history
    source_filename = fullfile(source_dir, jobdir, 'results.mat');
    if isfile(source_filename)
        fig = plotOptHistory(source_filename);
        % Save figure
        figbase = sprintf('history_%i', k);
        filename = fullfile(destin_dir, figbase);
        exportgraphics(fig, [filename, '.pdf'], 'Append', false);
        saveas(fig, [filename, '.png']);
        close(fig);
    end
end