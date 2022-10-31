%% Analyze results from buckling study 2
% Load job folders
D = dir('processed_data/buckling_study2/job*');

% Make struct of inputs iterating over all jobs
% Load first job to get fieldnames
dir1 = fullfile(D(1).folder, D(1).name);
inp1 = fullfile(dir1, 'input.mat');
dat1 = load(inp1); dat1.job = dir1;
fnames = fieldnames(dat1);
dcell = cell(size(fnames, 1));
% Making struct
input_data = cell2struct(dcell, fnames, 1);
input_data(1) = dat1;
% Iterating over all jobs
for k = 2:numel(D)
    dirk = fullfile(D(k).folder, D(k).name);
    inpk = fullfile(dirk, 'input.mat');
    datk = load(inpk);
    datk.job = dirk;
    input_data(k) = datk;
end
% Sort jobs by fieldnames
for k = 1:numel(fnames)-1
    % Check type of value and sort accordingly
    if isnumeric(input_data(1).(fnames{k}))
        [~, I] = sort([input_data.(fnames{k})], 'ascend');
    else
        [~, I] = sort({input_data.(fnames{k})});
    end
    % Permute input data
    input_data = input_data(I);
    D = D(I);
end
%% Reference data
Dref = dir('processed_data/buckling_study_reference1/job*');
fnames_ref = {...
    'Stats'; ...
    'xPhys'; ...
    'runtime'; ...
    'job'};
results_ref = cell2struct(cell(size(fnames_ref)), fnames_ref, 1);
for k = 1:numel(Dref)
    dirk = fullfile(Dref(k).folder, Dref(k).name);
    resk = load(fullfile(dirk, 'results.mat'), fnames_ref{1:end-1});

    resk.job = dirk;
    
    results_ref(k) = resk;
    fprintf('job(%03i): runtime %5.2f (min)\n', ...
        k, round(results_ref(k).runtime/60, 2));
end
%% Retrieve run data
fnames = {...
    'Stats'; ...
    'xPhys'; ...
    'runtime'; ...
    'speedup'; ...
    'job'};
results = cell2struct(cell(size(fnames)), fnames, 1);
for k = 1:numel(D)
    dirk = fullfile(D(k).folder, D(k).name);
    resk = load(fullfile(dirk, 'results.mat'), fnames{1:end-2});

    resk.speedup = results_ref(ceil(k/81)).runtime - resk.runtime;
    resk.job = dirk;
    
    results(k) = resk;
    fprintf('job(%03i): runtime %5.2f (min)\n', ...
        k, round(results(k).runtime/60, 2));
end

%% Analyze data
% Find maximum speedup
n0 = numel(D)/3;
speedups = [results.speedup];
max_speedups = zeros(3, 3);
for k = 1:3
    [S, I] = max(speedups((1+n0*(k-1)):(n0*k)));
    max_speedups(1:3, k) = [S, S/results_ref(k).runtime, I+n0*(k-1)];
end

% Plot correspodning designs next to reference design
ndesign = 2;
fig = figure();
for k = 1:3
    % load reference geometry
    load(fullfile(results_ref(k).job, 'results.mat'), ...
        'xPhys', 'i_img', 'j_img', 'sizex', 'sizey');
    % make figure of appropriate size
    dims = [sizex, sizey];
    [~, mindim] = min(dims);
    figsize = dims;
    figsize(mindim) = figsize(mindim)*ndesign;
    fig.Position(3:4) = figsize;
    ntiles = [1, 1];
    ntiles(mindim) = ntiles(mindim)+1;
    % make tiles
    t = tiledlayout(fig, ntiles(2), ntiles(1), ...
        'TileSpacing','tight');
    
    % plot reference design
    refimg = sparse(i_img, j_img, xPhys);
    nexttile(t, 1);
    imagesc(refimg);
    axis(gca, 'off');
    axis(gca, 'image')

    % plot rom design
    load(fullfile(results(max_speedups(3, k)).job, 'results.mat'), ...
        'xPhys')
    romimg = sparse(i_img, j_img, xPhys);
    nexttile(t, 2);
    imagesc(romimg);
    axis(gca, 'off');
    axis(gca, 'image')

    saveas(fig, sprintf('geometry%i.png', k));
end

%% Runtime of jobs compared to number of basis vectors
fig = figure();
t = tiledlayout(fig, 1, 3);

runs = 1:81+162;
params = {'NBASIS_STAT', 'NBASIS_EIGS', 'NBASIS_ADJT'};
for k = 1:3
    xdata = [input_data(runs).(params{k})];
    ydata = [results(runs).runtime];
    plot(nexttile(k), xdata, ydata, 'bo');
end

%% Set fixed number of basis for eigenvalue problem and same nb for stat and adjt
fig = figure();
