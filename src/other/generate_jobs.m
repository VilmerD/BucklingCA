function generate_jobs(jobDir, jobName, jobfarmParentdir, jobfarmDirnameBase, variables)
% GENERATE_JOBS(jobDir, jobName, jobfarmParentdir, jobfarmDirnameBase, variables) generates
% a jobfarm containing jobs with all combinations of the variables
% 

% Preprocess variables
variableNames = fieldnames(variables);
nVariables = numel(variableNames);
size_vec = ones(nVariables, 1);
for k = 1:nVariables
    size_vec(k) = numel(variables.(variableNames{k}));
end
nn = prod(size_vec);

% Make directory for tests. 
jobfarmDirname = jobfarmDirnameBase;
msg = zeros(1, 1);
num = 1;

% To avoid naming conflicts (ie creating a directory with the same name as an existing directory) 
% add a number at the end of the jobfarmDirname. Ex. if the directory jobfarm1 already exists create
% jobfarm2, or jobfarm3 until a unique name has been found
while ~isempty(msg)
    jobfarmDirname = sprintf([jobfarmDirnameBase, '%i'], num);
    [~, msg, ~] = mkdir(jobfarmParentdir, jobfarmDirname);
    num = num+1;
end

% Copy contents of scriptdir into the jobfarm directory
copyfile(jobDir, fullfile(jobfarmParentdir, jobfarmDirname));

% Make input files
jobvariables = struct();
subs = cell(1, nVariables);
for k = 1:nn
    % Make job dir
    jobDirname = sprintf('job_%i', k-1);
    fullJobDirname = fullfile(jobfarmParentdir, jobfarmDirname, jobDirname);
    mkdir(fullJobDirname);
    
    % Copy job to jobdir
    copyfile(fullfile(jobfarmParentdir, jobfarmDirname, jobName), ...
        fullfile(fullJobDirname, 'job.m'));
    
    % Make input data
    [subs{:}] = ind2sub(size_vec, k);
    for i = 1:numel(variableNames)
        varnamei = variableNames{i};
        vari = variables.(varnamei){subs{i}};
        jobvariables.(varnamei) = vari;
    end
    
    % Save input data
    dataname = fullfile(fullJobDirname, 'input.mat');
    save(dataname, '-struct', 'jobvariables');
end

end