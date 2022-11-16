%% Evaluate work done in terms of linear solves and factorizations
% Load dummy data
resfile = 'processed_data/rom13/job_0/results.mat';
load(resfile, 'profile_data');
fnc_table = profile_data.FunctionTable;
%% How many triangular solves in eigs?
% How many triangular solves were excecuted when solving eigenproblem
trisolver = 'eigs>@(v)solve(dR,permB''*applyA(permB*solve(dR,v,false)),true)';  % name of function
itrisolver = cellfun(@(n) strcmp(n, trisolver), {fnc_table.FunctionName});      % index in table
ncalls_eigs = fnc_table(itrisolver).NumCalls;                                   % number of calls

% Triangluar solves in CAEEON
name_caeeon = 'CAEEON';
lines_caeeon = [12, 24];
tri_caeeon = countCalls(fnc_table, name_caeeon, lines_caeeon)

%% How many triangular solves for static problem / adjoints?
name_msolveq = 'msolveq';
lines_solveq = 49;
tri_solveq = countCalls(fnc_table, name_msolveq, lines_solveq)

name_casbon = 'CASBON';
lines_casbon = [9, 18];
tri_casbon = countCalls(fnc_table, name_casbon, lines_casbon)

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
