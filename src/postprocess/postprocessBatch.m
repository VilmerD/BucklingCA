%% Given a design computes blf, volume fraction, compliance etc
farm = 'ref010';

source_dir = fullfile('processed_data/', farm);

datmat_des = struct('VF', {}, 'C', {}, 'L', {}, 'pL', {});
datmat_wrk = struct('nIts', {}, 'nFact', {}, 'nCA', {}, 'nTri', {}, ...
    'tsol', {}, 'wT', {}, 'wTT', {});

nj = numel(dir(fullfile(source_dir, 'job_*')));
for k = 0:nj-1
    jobdir = fullfile(source_dir, sprintf('job_%i', k));
    if isfile(fullfile(jobdir, 'results.mat'))
        datmat_des(k+1) = evaluateDesign(jobdir);
        datmat_wrk(k+1) = evaluateOpt(jobdir);
    end

    if ~mod(k, nj/10); fprintf('%%'); end
end
fprintf("\n");

%% Save data
destindir = 'processed_data/processed_batches/';
datfilename = fullfile(destindir, sprintf('%s.mat', farm));
S = struct(sprintf('design_%s', farm), datmat_des, sprintf('work_%s', farm), datmat_wrk);
save(datfilename, '-struct', 'S');
