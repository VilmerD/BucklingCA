%% Plotting designs of reference problem using exact analysis
jobs = dir('processed_data/buckling_study_reference1/job*');

nameform = 'results/images/desgins/%s_reference.png';
for k = 1:3
    results_file = fullfile(jobs(k).folder, jobs(k).name, 'results.mat');
    load(results_file, 'xPhys', 'i_img', 'j_img', 'domain');
    outfile = sprintf(nameform, domain);
    plotDesignAndSave(xPhys, i_img, j_img, outfile);
end