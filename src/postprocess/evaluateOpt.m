function S = evaluateOpt(inpfile, resfile)
% Counts number of factorizations etc for solution
load(resfile, 'Stats');

% Iterations and factorizations
nIts    = sum(Stats(:, 1)  ~=0);
nFact   = sum(Stats(:, end)==1);
nCA     = nIts - nFact;

% Triangular solves
load(resfile, 'NBASIS_STAT', 'NBASIS_EIGS', 'NBASIS_ADJT', 'lambda');
nevals = numel(lambda);
tpi_stat_fc = 1;
tpi_eigs_fc = nevals;
tpi_adjt_fc = nevals;
tpi_stat_ca = NBASIS_STAT;
tpi_eigs_ca = NBASIS_EIGS*2 + 2*(nevals-2);
tpi_adjt_ca = NBASIS_ADJT*2 + 2*(nevals-2);
if isnan(NBASIS_STAT) || isnan(NBASIS_EIGS) || isnan(NBASIS_ADJT)
    nTri = [...
        nFact*tpi_stat_fc...
        nFact*tpi_eigs_fc...
        nFact*tpi_adjt_fc];
else
    nTri = [...
        nCA*tpi_stat_ca + nFact*tpi_stat_fc...
        nCA*tpi_eigs_ca + nFact*tpi_eigs_fc...
        nCA*tpi_adjt_ca + nFact*tpi_adjt_fc];
end
nTri(4) = sum(nTri);

% Wall - time
load(resfile, 'runtime', 'tsol');
wT = sum(tsol, 1);
% Save in struct
datnames = {'nIts', 'nFact', 'nCA', 'nTri', 'tsol', 'wT', 'wTT', 'runtime'};
datvals = {nIts, nFact, nCA, nTri, tsol, wT, sum(wT), runtime};
S = cell2struct(datvals, datnames, 2);
end