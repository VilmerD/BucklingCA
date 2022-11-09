function S = evaluateDesign(inpfile, resfile)
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

% Save in struct
datnames = {'VF', 'C', 'L', 'pL', 'Ls'};
datvals = {VF, C, L, pL, Ls};
S = cell2struct(datvals, datnames, 2);
end