%% Commit Jobs
% This script moves selected jobs to report directory and collects work and
% design data to be sent to LaTeX document. 

%% Select data
srcdirs = {'designs', 'history'};
dstdir = 'report';

% Batches and jobs to commit
farms = {'ref022', 'rom034'};
jobs = {0:2, 9:11};

% Load safety-factors
geoms = cell(numel(farms), 1);
sfs = zeros(numel(farms), 1);
for k = 1:numel(farms)
    fk = farms{k};
    jk = jobs{k};
    for l = 1:numel(jk)
        load(fullfile(...
            'processed_data/', ...
            fk, sprintf('job_%i', jk(l)), ...
            'input.mat'), ...
            'sf', 'domain');
        sfs(k, l) = sf;
        geoms{k, l} = domain;
    end
end

%% Format table
% Format LaTeX table from data, because I'm too lazy to do data entry
filename = 'results/report/values.txt';
fid = fopen(filename, 'w+');

% Write jobs and numbers in jobs
fprintf(fid, 'JOBS\n\n');
fprintf(fid, '%-10s \t %-10s\n', 'Farm names', 'Job IDs');
for k = 1:length(jobs)
    fk = farms{k};
    jk = jobs{k};
    format = ['%-10s \t ', repmat('%i, ', 1, numel(jk)), '\n'];
    fprintf(fid, format, fk, jk);
end
fprintf(fid, '\n');

%% Move figures
% Pattern to remove numbering
pat = '(?<case>\D+)\d+(?<ext>\D+)';

% Loop over farms (batches)
for k = 1:numel(farms)
    fk = farms{k};
    jk = jobs{k};

    % Loop over designs and history
    for p = 1:numel(srcdirs)
        % Make dir
        srck = srcdirs{p};
        srcdir = fullfile('results', srck);
        dk = fullfile('results', dstdir, srck, fk);
        if isfolder(dk); rmdir(dk, 's'); end
        mkdir(dk);

        % Loop over jobs
        for j = 1:3
            src = fullfile('results', srck, ...
                sprintf('blfopt_%s', fk), ...
                sprintf('*%i*', jk(j)));
            dst = dk;
            copyfile(src, dst);
        end

        % Remove numbers (for LaTeX format purposes)
        ck = dir(dk);
        for j = (3:numel(ck))
            fj = ck(j).name;
            names = regexp(fj, pat, 'names');
            newname = [names.case, names.ext];
            movefile(fullfile(ck(j).folder, ck(j).name), ...
                fullfile(ck(j).folder, newname));
        end
    end
end

%% Compliance reference
lsdat = zeros(size(sfs));
for k = 1:numel(farms)
    jk = jobs{k};
    for l = 1:numel(jk)
        load(fullfile('data/compliance_reference', geoms{k, l}), 'lambda');
        lsdat(k, l) = lambda(1)*sfs(k);
    end
end

%% Design data
fprintf(fid, 'DESIGN DATA\n\n');
desdat = cell(2, 1);
for k = 1:numel(farms)
    S = load(fullfile('processed_data/processed_batches', farms{k}));
    desdatall = S.(sprintf('design_%s', farms{k}));
    desdat{k} = desdatall(jobs{k}+1);
end

mform = '%-8s & %8.4f & %8.4f & %8.4f & %8.4f \\\\\n';
for k = 1:numel(jobs{1})
    for l = 1:numel(farms)
        datj = desdat{l};
        datjk = datj(k);
        fj = farms{l};
        fprintf(fid, mform, upper(fj(1:3)), datjk.VF, datjk.C/1e3, ...
            lsdat(l, k)/datjk.L(1), lsdat(l, k)/datjk.pL);
    end
end

fprintf(fid, '\n');
%% Work data
fprintf(fid, 'WORK DATA\n\n');
wrkdat = cell(2, 1);
for k = 1:numel(farms)
    S = load(fullfile('processed_data/processed_batches', farms{k}));
    workdatall = S.(sprintf('work_%s', farms{k}));
    wrkdat{k} = workdatall(jobs{k}+1);
end

mform = '%-8s & %8i & %8i & %8i & %8i \\\\ \n';
for k = 1:numel(jobs{1})
    for j = 1:numel(farms)
        datj = wrkdat{j};
        datjk = datj(k);
        fj = farms{j};
        fprintf(fid, mform, upper(fj(1:3)), datjk.nFact, datjk.nTri(1), datjk.nTri(2), datjk.nTri(3));
    end
end

%% Finish
fclose(fid);