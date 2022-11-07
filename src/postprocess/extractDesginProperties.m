%% Given a design computes blf, volume fraction, compliance etc
% Assuming square elements (ie L = H)
% Define data (temporary)
farm = 'rom3';
job = 'job_8';
datfile_inp = fullfile('processed_data', farm, job, 'input.mat');
datfile_res = fullfile('processed_data', farm, job, 'results.mat');
% Load data
load(datfile_inp, 'domain');
load(datfile_res, 'xPhys', ...
    'Emin', 'Emax', 'nu',  'pE', 'pS', 'KE');
switch lower(domain)
    case 'column'
        sizex = 240;
        sizey = 120;
    case 'spire'
        sizex = 240;
        sizey = 120;
    case 'twobar'
        sizex = 120;
        sizey = 280;
end
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
vf = sum(xPhys)*helem^2/(sizex*sizey);

% Compute system matrices
iK = reshape(kron(edofMat, ones(8,1))', 64*nelem, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelem, 1);
% Compute penalized density field
xPenL = Emin + (Emax - Emin)*xPhys.^pE;
xPenG = Emax*xPhys.^pE;

% Compute linear first
sK = kron(xPenL, KE(:));
K = sparse(iK, jK, sK);
K = (K' + K)/2;
[R, CHOLOK, P] = chol(K(freedofs, freedofs), 'matrix');
U = msolveq(K, F, bc, R, P);

% Then compute nonlinear
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       % constitutive matrix - plane stress
sG = computeKG(U, D, edofMat, helem, xPenG);
KNL = sparse(iK, jK, sG);
KNL = (KNL' + KNL)/2;

% Compute comp
comp = U'*F;

% Compute blf
[evecs, evals] = meigenSM(-KNL, K, bc, 1, R, P);
blf = 1/evals(1);

% Display
fprintf('VQ  : %+10.3f \nCOMP: %10.3f \nBLF : %+10.3f\n', vf, comp, blf);