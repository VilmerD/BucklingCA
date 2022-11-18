function S = evaluateDesign(resfile)
% Load data
load(resfile, 'INPVARS', 'stats', 'xPhys');
% Load domain data
load(fullfile('data/domain_data', sprintf('%s.mat', INPVARS.domain)), ...
    'sizex', 'sizey');
% Generate Geometry
genfun = str2func(sprintf('generate_%s', INPVARS.domain));
[X,T,~,~,~,~,F,freedofs] = genfun(sizex,sizey,INPVARS.helem,0);
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
VF = sum(xPhys)*INPVARS.helem^2/(sizex*sizey);
% Compute system matrices
iK = reshape(kron(edofMat, ones(8,1))', 64*nelem, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelem, 1);
% Compute penalized density field
Emax = 2e5;
Emin = Emax*1e-6;
pE = stats(end, 14);
pS = pE;
xPenL = Emin + (Emax - Emin)*xPhys.^pE;
xPenG = Emax*xPhys.^pS;
% Compute linear first
nu = 0.3;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
sK = kron(xPenL, KE(:));
K = sparse(iK, jK, sK);
K = (K' + K)/2;
[R, ~, P] = chol(K(freedofs, freedofs), 'matrix');
U = msolveq(K, F, bc, R, P);
% Then compute nonlinear
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       % constitutive matrix - plane stress
sG = computeKG(U, D, edofMat, INPVARS.helem, xPenG);
KNL = sparse(iK, jK, sG);
KNL = (KNL' + KNL)/2;
% Compute comp
C = U'*F;
% Compute blf
pN = stats(end, 15);
[~, evals] = meigenSM(-KNL, K, bc, INPVARS.nevals, R, P);
L = 1./diag(evals);
pL = 1./norm(1./L, pN);
% Target blf
load(fullfile('data/compliance_reference/', sprintf('%s.mat', INPVARS.domain)), ...
    'lambda');
Ls = INPVARS.sf*lambda(1);
% Display
fprintf('VQ  : %+10.3f \nCOMP: %10.3f \nBLF : %+10.3f\n', VF, C, L(1));
% Save in struct
datnames = {'VF', 'C', 'L', 'pL', 'Ls'};
datvals = {VF, C, L, pL, Ls};
S = cell2struct(datvals, datnames, 2);
end