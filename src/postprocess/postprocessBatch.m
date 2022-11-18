%% Given a design computes blf, volume fraction, compliance etc
% Assuming square elements (ie L = H)
% Define data (temporary)
farm = 'rom15';
jobnums = 0:7;
%% Extract design data
datmat_des = struct('VF', {}, 'C', {}, 'L', {}, 'pL', {}, 'Ls', {});
datmat_wrk = struct('nIts', {}, 'nFact', {}, 'nCA', {}, 'nTri', {}, ...
    'tsol', {}, 'wT', {}, 'wTT', {});
for k = 1:numel(jobnums)
    jobk = sprintf('job_%i', jobnums(k));
    datfile_res = fullfile('processed_data', farm, jobk, 'results.mat');
    if isfile(datfile_res)
        datmat_des(k) = evaluateDesign(datfile_res);
        datmat_wrk(k) = evaluateOpt(datfile_res);
    end
end

%% Save data
destindir = 'processed_data/processed_batches/';
datfilename = fullfile(destindir, sprintf('%s.mat', farm));
S = struct(sprintf('design_%s', farm), datmat_des, sprintf('work_%s', farm), datmat_wrk);
save(datfilename, '-struct', 'S');
