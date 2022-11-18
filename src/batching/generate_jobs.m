function generate_jobs(jobDir, jobName, jobfarmParentdir, jobfarmDirnameBase, variables)
% GENERATE_JOBS(jobDir, jobName, jobfarmParentdir, jobfarmDirnameBase, variables) generates
% a jobfarm containing jobs with all combinations of the variables
%

% Preprocess variables
variables_permuted = permuteVariables(variables);

% Make directory for tests.
jobfarmDirname = jobfarmDirnameBase;
msg = zeros(1, 1);
num = 1;

% To avoid naming conflicts (ie creating a directory with the same name as an existing directory)
% add a number at the end of the jobfarmDirname. Ex. if the directory jobfarm1 already exists create
% jobfarm2, or jobfarm3 until a unique name has been found
while ~isempty(msg)
    jobfarmDirname = sprintf([jobfarmDirnameBase, '%03i'], num);
    [~, msg, ~] = mkdir(jobfarmParentdir, jobfarmDirname);
    num = num+1;
end

% Copy contents of scriptdir into the jobfarm directory
copyfile(fullfile(jobDir, jobName), fullfile(jobfarmParentdir, jobfarmDirname));
copyfile(fullfile(jobDir, 'master_script.sh'), fullfile(jobfarmParentdir, jobfarmDirname));
copyfile(fullfile(jobDir, 'worker_script.sh'), fullfile(jobfarmParentdir, jobfarmDirname));

% Make input files
for k = 1:numel(variables_permuted)
    % Make job dir
    jobDirname = sprintf('job_%i', k-1);
    fullJobDirname = fullfile(jobfarmParentdir, jobfarmDirname, jobDirname);
    mkdir(fullJobDirname);

    % Copy job to jobdir
    copyfile(fullfile(jobfarmParentdir, jobfarmDirname, jobName), ...
        fullfile(fullJobDirname, 'job.m'));

    % Save input data
    inputk = variables_permuted(k);
    dataname = fullfile(fullJobDirname, 'input.mat');
    save(dataname, '-struct', 'inputk');
end

end

function variables_permuted = permuteVariables(variables)
% Get number of permutations 
variableNames = fieldnames(variables);
nVariables = numel(variableNames);
size_vec = ones(nVariables, 1);
for k = 1:nVariables
    size_vec(k) = numel(variables.(variableNames{k}));
end
nn = prod(size_vec);

% Make a struct for each permutation
jobvariables = struct();
subs = cell(1, nVariables);
variables_permuted = cell2struct(...
    cell(size(variableNames)), variableNames, find(max(size(variableNames))));
for k = 1:nn
    [subs{:}] = ind2sub(size_vec, k);
    for i = 1:numel(variableNames)
        varnamei = variableNames{i};
        vari = variables.(varnamei){subs{i}};
        jobvariables.(varnamei) = vari;
    end
    variables_permuted(k) = jobvariables;
end
end