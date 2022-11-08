%% Given a design computes blf, volume fraction, compliance etc
% Assuming square elements (ie L = H)
% Define data (temporary)
farm = 'ref7';
jobnums = 0:11;
%% Extract design data
if isempty(gcp); pool = parpool(6); end
datmat_des = struct(...
    'VF', {}, ...
    'C', {}, ...
    'L', {}, ...
    'pL', {}, ...
    'Ls', {});
datmat_wrk = struct(...
    'nIts', {}, ...
    'nFact', {}, ...
    'nCA', {}, ...
    'nTri', {}, ...
    'tsol', {}, ...
    'wT', {}, ...
    'wTT', {}, ...
    'runtime', {});
parfor (k = 1:numel(jobnums), pool)
    jobk = sprintf('job_%i', jobnums(k));
    datfile_inp = fullfile('processed_data/batches', farm, jobk, 'input.mat');
    datfile_res = fullfile('processed_data/batches', farm, jobk, 'results.mat');
    [VFk, Ck, Lk, pLk, Lsk] = extractDesignProperties(datfile_inp, datfile_res);
    datmat_des(k) = struct('VF', VFk, 'C', Ck, 'L', Lk, 'pL', pLk, 'Ls', Lsk);
    datmat_wrk(k) = extractWork(datfile_inp, datfile_res);
end

%% Save data
destindir = 'processed_data/processed_batches/';
datfilename = fullfile(destindir, sprintf('%s.mat', farm));
save(datfilename, 'datmat_des', 'datmat_wrk');

%%
function [VF, C, L, pL, Ls] = extractDesignProperties(inpfile, resfile)
% Load data
load(inpfile, 'domain');
load(resfile, 'sizex', 'sizey', ...
    'Emin', 'Emax', 'nu',  'pE', 'pS', 'KE', ...
    'xPhys');
helem = sqrt(sizex*sizey/numel(xPhys));
% Generate Geometry
genfun = str2func(sprintf('generate_%s', domain));
[X,T,i_img,j_img,soldis,voids,F,freedofs] = genfun(sizex,sizey,helem,0);
nelem = size(T, 1);
ndof = 2*size(X,1);
edofMat = zeros(nelem,8);
edofMat(:,1:2) = [2*T(:,1)-1 2*T(:,1)];
edofMat(:,3:4) = [2*T(:,2)-1 2*T(:,2)];
edofMat(:,5:6) = [2*T(:,3)-1 2*T(:,3)];
edofMat(:,7:8) = [2*T(:,4)-1 2*T(:,4)];
% BC
bc = (1:ndof)';
bc(freedofs) = [];
bc(:, 2) = 0;

% Compute volume fraction
VF = sum(xPhys)*helem^2/(sizex*sizey);

% Compute system matrices
iK = reshape(kron(edofMat, ones(8,1))', 64*nelem, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelem, 1);
% Compute penalized density field
xPenL = Emin + (Emax - Emin)*xPhys.^pE;
xPenG = Emax*xPhys.^pS;

% Compute linear first
sK = kron(xPenL, KE(:));
K = sparse(iK, jK, sK);
K = (K' + K)/2;
[R, ~, P] = chol(K(freedofs, freedofs), 'matrix');
U = msolveq(K, F, bc, R, P);

% Then compute nonlinear
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       % constitutive matrix - plane stress
sG = computeKG(U, D, edofMat, helem, xPenG);
KNL = sparse(iK, jK, sG);
KNL = (KNL' + KNL)/2;

% Compute comp
C = U'*F;

% Compute blf
load(resfile, 'lambda', 'pN')
nevals = numel(lambda);
[~, evals] = meigenSM(-KNL, K, bc, nevals, R, P);
L = 1./diag(evals);
pL = 1./norm(1./L, pN);

% Target blf
load(resfile, 'lamstar');
Ls = lamstar;

% Display
fprintf('VQ  : %+10.3f \nCOMP: %10.3f \nBLF : %+10.3f\n', VF, C, L(1));
end

function S = extractWork(inpfile, resfile)
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


% Wall - time
load(resfile, 'runtime', 'tsol');
wT = sum(tsol, 1);
% Save in struct
datnames = {'nIts', 'nFact', 'nCA', 'nTri', 'tsol', 'wT', 'wTT', 'runtime'};
datvals = {nIts, nFact, nCA, nTri, tsol, wT, sum(wT), runtime};
S = cell2struct(datvals, datnames, 2);
end