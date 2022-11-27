function S = evaluateOpt(jobdir)
% Counts number of factorizations etc for solution
load(fullfile(jobdir, 'results.mat'), 'stats', 'profile_data', 'tsol');

% Iterations and factorizations
nIts    = sum(stats(:, 1)  ~=0);
nFact   = sum(stats(:, end)==1);
nCA     = nIts - nFact;

fnc_table = profile_data.FunctionTable;

%%% Triangular solves %%%
% lines are the lines in the source code where mldivide is called, this
% must be done since matlab does not log mldivide on its own

% for linear equilibrium and adjoints
nTri(1) = countCalls(fnc_table, 'msolveq', 49);
nTri(1) = nTri(1) + countCalls(fnc_table, 'CASBON', [11, 19]);

% for eigenvalueproblem
nTri(2) = countCalls(fnc_table, 'CAEEON', [11, 22]);
% the triangular solves in eigs is in a lambda-expression
trisolver = 'eigs>@(v)solve(dR,permB''*applyA(permB*solve(dR,v,false)),true)';  % name of function
itrisolver = cellfun(@(n) strcmp(n, trisolver), {fnc_table.FunctionName});      % index in table
nTri(2) = nTri(2) + fnc_table(itrisolver).NumCalls;         

% total number of trisolves
nTri(end+1) = sum(nTri);

% walltime
wT = sum(tsol, 1);

% Save data
datnames = {'nIts', 'nFact', 'nCA', 'nTri', 'tsol', 'wT', 'wTT'};
datvals = {nIts, nFact, nCA, nTri, tsol, wT, sum(wT)};
S = cell2struct(datvals, datnames, 2);
end

function ncalls = countCalls(fnc_table, fnc_name, lines)
indx = cellfun(@(n) strcmp(n, fnc_name), {fnc_table.FunctionName});
fnc = fnc_table(indx);
ncalls = 0;
if ~isempty(fnc)
    exelines = fnc.ExecutedLines;
    for l = lines
        ncalls = ncalls + exelines(exelines(:, 1)==l, 2);
    end
end
end