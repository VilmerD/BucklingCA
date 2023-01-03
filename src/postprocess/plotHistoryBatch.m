%% Plot the compliance and blf's
farm = 'rom034';

source_dir = fullfile('processed_data', farm);
destin_dir = fullfile('results/history', sprintf('blfopt_%s', farm));
mkdir(destin_dir);

nj = numel(dir(fullfile(source_dir, 'job_*')));
for k = 0:nj-1
    jobdir = sprintf('job_%i', k);
    % Plot history
    source_filename = fullfile(source_dir, jobdir, 'results.mat');
    if isfile(source_filename)
        % Make figure
        fig = plotOptHistory(source_filename);
        % Save figure
        load(source_filename, 'INPVARS');
        filename = fullfile(destin_dir, sprintf('%s%i', INPVARS.domain, k));
        exportgraphics(fig, [filename, '.pdf'], 'Append', false);
        saveas(fig, [filename, '.png']);
        close(fig);
    end
end