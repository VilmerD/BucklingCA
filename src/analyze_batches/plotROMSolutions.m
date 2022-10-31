%% Plotting designs of reference problem using exact analysis
jobs = dir('processed_data/buckling_study2/job*');
params = struct('NBASIS_STAT', 8, 'NBASIS_EIGS', 2, 'NBASIS_ADJT', 8, ...
    'CHOL_UPDATE_PACE', 10);
fn = fieldnames(params);

nameform = 'results/images/designs/%s_ca.png';
for k = 1:numel(jobs)
    input_file = fullfile(jobs(k).folder, jobs(k).name, 'input.mat');
    inpts = load(input_file, fn{:});
    
    % check if inputs correspond to params
    plot_design = 1;
    for l = 1:numel(fn)
        if inpts.(fn{l}) ~= params.(fn{l}); plot_design = 0; break; end
    end

    % read result data
    if plot_design
        results_file = fullfile(jobs(k).folder, jobs(k).name, 'results.mat');
        load(results_file, 'xPhys', 'i_img', 'j_img', 'domain');
        outfile = sprintf(nameform, domain);
        plotDesignAndSave(xPhys, i_img, j_img, outfile);
    end
    k
end