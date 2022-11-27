% Script used to format LaTeX table from data, because I'm too lazy to do 
% it myself so I make MATLAB do it instead :)

%% Load data
file = 'rom014';
indx = 55:57;
S = load(fullfile('processed_data/processed_batches', file));
%% design
designdat_all = S.(sprintf('design_%s', file));
designdat = designdat_all(indx);

mform = '& %5.4f & %5.4f & %5.4f & %5.4f & \n';
for k = 1:numel(designdat)
    datk = designdat(k);
    fprintf(mform, datk.VF, datk.C/1e4, datk.L(1), datk.pL);
end

%% work
workdat_all = S.(sprintf('work_%s', file));
workdat = workdat_all(indx);

mform = '& %i & %i & %i & %i \\\\ \n';
for k = 1:numel(workdat)
    datk = workdat(k);
    fprintf(mform, datk.nFact, datk.nTri(1), datk.nTri(2), datk.nTri(3));
end