%% Given a design computes blf, volume fraction, compliance etc
% Assuming square elements (ie L = H)
% Define data (temporary)
farm = 'ref8';
jobnums = 0:8;
%% Extract design data
pool = gcp;
if isempty(pool); pool = parpool(6); end
datmat_des = struct('VF', {}, 'C', {}, 'L', {}, 'pL', {}, 'Ls', {});
datmat_wrk = struct('nIts', {}, 'nFact', {}, 'nCA', {}, 'nTri', {}, ...
    'tsol', {}, 'wT', {}, 'wTT', {}, 'runtime', {});
parfor (k = 1:numel(jobnums), pool)
    jobk = sprintf('job_%i', jobnums(k));
    datfile_inp = fullfile('processed_data/batches', farm, jobk, 'input.mat');
    datfile_res = fullfile('processed_data/batches', farm, jobk, 'results.mat');
    datmat_des(k) = evaluateDesign(datfile_inp, datfile_res);
    datmat_wrk(k) = evaluateOpt(datfile_inp, datfile_res);
end

%% Save data
destindir = 'processed_data/processed_batches/';
datfilename = fullfile(destindir, sprintf('%s.mat', farm));
S = struct(sprintf('design_%s', farm), datmat_des, sprintf('work_%s', farm), datmat_wrk);
save(datfilename, '-struct', 'S');
