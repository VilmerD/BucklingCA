% Define job
jobdir = 'processed_data/ref012/job_0';

% Get files
INPVARS = load(fullfile(jobdir, 'input.mat'));
RESVARS = load(fullfile(jobdir, 'results.mat'));
% Load domain data
helem = INPVARS.helem;
load(fullfile('data/domain_data', sprintf('%s.mat', INPVARS.domain)), ...
    'sizex', 'sizey');
% Generate Geometry
genfun = str2func(sprintf('generate_%s', INPVARS.domain));
[X,T,~,~,~,~,F,freedofs] = genfun(sizex,sizey,helem,0);
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
xPhys = RESVARS.xPhys;
VF = sum(xPhys)*helem^2/(sizex*sizey);
% Compute system matrices
iK = reshape(kron(edofMat, ones(8,1))', 64*nelem, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelem, 1);
% Compute penalized density field
Emin = INPVARS.mpara(1); Emax = INPVARS.mpara(2);
pE = INPVARS.pphysend; pS = pE;
xPenL = Emin + (Emax - Emin)*xPhys.^pE;
xPenG = Emax*xPhys.^pS;
%% Compute linear stiffness matrix
nu = INPVARS.mpara(3);
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
sK = kron(xPenL, KE(:));
K = sparse(iK, jK, sK); K = (K' + K)/2;
% Solve linear equilibrium
[R, ~, P] = chol(K(freedofs, freedofs), 'matrix');
U = msolveq(K, F, bc, R, P);
%% Compute nonlinear stiffness matrix
% constitutive matrix - plane stress
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       
% Domain data
% Base matrices
B = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    0 -1/2 0 -1/2 0 1/2 0 1/2
    -1/2 -1/2 -1/2 1/2 1/2 1/2 1/2 -1/2];           % strain-displacement matrix
BNL = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    -1/2 0 -1/2 0 1/2 0 1/2 0
    0 -1/2 0 1/2 0 1/2 0 -1/2
    0 -1/2 0 -1/2 0 1/2 0 1/2];                     % BNL matrix
% Compute stress and strain
EPS = U(edofMat)*B';
SIG = EPS*D;
% Compute element matrices
sG = zeros(nelem*64,1);
for el = 1:nelem
    tau = [SIG(el,1) SIG(el,3) 0 0;
            SIG(el,3) SIG(el,2) 0 0 ;
            0 0 SIG(el,1) SIG(el,3);
            0 0 SIG(el,3) SIG(el,2)];
        BNLTtau = BNL'*tau;
        BNLTtauBNL = BNLTtau*BNL;
        l1 = (el-1)*64+1; l2 = el*64;
        sG(l1:l2) = xPenG(el)*helem^2*BNLTtauBNL(:);
end
KNL = sparse(iK, jK, sG); KNL = (KNL' + KNL)/2;
% Compute blf
[evecs, evals] = meigenSM(-KNL, K, bc, INPVARS.nevals, R, P);
evecs = evecs./sqrt(dot(evecs, K*evecs));

%% Compute nodal positions
enodMat = (edofMat(:, 1:2:end) + 1)/2;
XC = X(:, 1); 
YC = X(:, 2);
ex = XC(enodMat);
ey = YC(enodMat);

%% Compute nodal displacements
k = 1;
phik = evecs(:, k);
ed = phik(edofMat);

exd = ex+
