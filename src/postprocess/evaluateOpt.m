function S = evaluateOpt(resfile)
% Counts number of factorizations etc for solution
load(resfile, 'stats');

% Iterations and factorizations
nIts    = sum(stats(:, 1)  ~=0);
nFact   = sum(stats(:, end)==1);
nCA     = nIts - nFact;

load(resfile, 'profile_data');
fnc_table = profile_data.FunctionTable;

% Triangular solves
name_msolveq = 'msolveq';
lines_solveq = 49;
nTri(1) = countCalls(fnc_table, name_msolveq, lines_solveq);

name_casbon = 'CASBON';
lines_casbon = [9, 18];
nTri(1) = nTri(1) + countCalls(fnc_table, name_casbon, lines_casbon);

name_caeeon = 'CAEEON';
lines_caeeon = [12, 24];
nTri(2) = countCalls(fnc_table, name_caeeon, lines_caeeon);

trisolver = 'eigs>@(v)solve(dR,permB''*applyA(permB*solve(dR,v,false)),true)';  % name of function
itrisolver = cellfun(@(n) strcmp(n, trisolver), {fnc_table.FunctionName});      % index in table
nTri(2) = nTri(2) + fnc_table(itrisolver).NumCalls;         

nTri(end+1) = sum(nTri);

% Wall - time
load(resfile, 'tsol');
wT = sum(tsol, 1);

% Save in struct
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