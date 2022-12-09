%% Copy compliance reference data to destination
bdir = 'compref006';
source_dir = fullfile('processed_data/', bdir);
destin_dir = 'data/compliance_reference';
for k = 0:2
    jobdir = fullfile(source_dir, sprintf('job_%i', k));
    inpfile = fullfile(jobdir, 'input.mat');
    load(inpfile, 'domain');
    dest = [domain, '.mat'];
    resfile = fullfile(jobdir, 'results.mat');
    copyfile(resfile, fullfile(destin_dir, dest))
end