%% Minimum max mu (1/lambda) s.t. volume
% Approximation of max mu using p-norm
% Single density field with eta=0.5
%% Add neccessary paths
addpath(genpath('~/Documents/Projects/BucklingCA/data/'));
addpath(genpath('~/Documents/Projects/BucklingCA/src/generate_geometry'));
addpath(genpath('~/Documents/Projects/BucklingCA/src/solvers'));
%% Display some information to user
fprintf('Starting time of job: %s\n', datestr(datetime('now')));
fprintf('Working directory: %s\n', pwd);
tstart = tic;
profile on;
%% Load run data
% If there is an inputfile named input.mat then load data from it, otherwise
% use default parameters
if isfile(fullfile(cd, 'input.mat'))
    variables = 'input.mat';
else
    variables = 'data/default_data.mat';
end
fprintf('Input variables: \n');
details(load(variables));
%% Domain size and discretization
% element size (all elements are square)
load(variables, 'domain', 'helem');
%% Prepare design domain
% X = nodal coordinates
% T = element connectivity and data
% i_img,j_img = for displaying design as a matrix using imagesc
% Load domain size and volume fraction
domaindata = fullfile('data/domain_data', sprintf('%s.mat', domain));
load(domaindata, 'sizex', 'sizey', 'volfrac');
% Generate geometry
geometry_generation = str2func(sprintf('generate_%s', domain));
[X,T,i_img,j_img,solids,voids,F,freedofs] = geometry_generation(sizex,sizey,helem,0);
%% Buckling parameters
load(variables, 'nevals');          % number of eigenvalues to consider
compref = fullfile('data/compliance_reference', sprintf('%s.mat', domain));
load(compref, 'lambda');
lamstar = lambda(1);                % Critical BLF (Lc)
load(variables, 'sf');              % Safety factor                
lamstar = lamstar*sf;               % Target BLF is safety factor * Lc
mustar = 1/lamstar;
%% Continuation control
load(variables, 'minloop', 'maxloop')

load(variables, 'pphysstart', 'pphysend', 'dpphys')
contsteps = ceil((pphysend-pphysstart)/dpphys);
pE = pphysstart;
pS = pphysstart;

load(variables, 'pNstart', 'pNend');
dpN = (pNend - pNstart)/contsteps;
pN = pNstart;

load(variables, 'betastart', 'betaend');
dbeta = (betaend - betastart)/contsteps;
beta = betastart;

load(variables, 'contstart')
pace = (minloop-contstart)/contsteps;
contnext = contstart;

loop = 0;
%% Track solution accuracy and work
stats = zeros(minloop,11+2+nevals);
tsol = zeros(maxloop, 4);
res_stat = zeros(minloop, 1);
res_eigs = zeros(minloop, nevals);
res_adjt = zeros(minloop, nevals);
mag_stat = res_stat;
mag_eigs = res_eigs;
mag_adjt = res_adjt;
%% Initialize CA
% Booleans controlling whether or not to use CA for
% the individual problems
SOLVE_STAT_EXACTLY = 0;
SOLVE_EIGS_EXACTLY = 0;
SOLVE_ADJT_EXACTLY = 0;
% Number of basis vectors
load(variables, 'NBASIS_STAT', 'NBASIS_EIGS', 'NBASIS_ADJT');
% Orthogonalization type for eigenproblem
load(variables, 'CAorthotype');
CAorthovecs = [];          % Old eigenmodes if orthotype is 'old', else empty
% When CA should be started, and how often to factorize
load(variables, 'CA_START', 'CHOL_UPDATE_PACE', 'CA_ANG_TOL')
%% Initialize figure
nc = 1;
nr = 3;
load(variables, 'SHOW_DESIGN')
if SHOW_DESIGN
    fig = figure(...
        'MenuBar', 'none');
    for c = 1:nc
        for r = nr:-1:1
            i = (c-1)*nr + r;
            ax(i) = axes(fig, ...
                'Units', 'Normalized', ...
                'Position', [(c-1)/nc, (nr-r)/nr, 1/nc, 1/nr]);
            hold(ax(i), 'on');
            axis(ax(i), 'off');
            axis(ax(i), 'tight');
            colormap(ax(i), 'summer');
            pause(0.05);
        end
    end
    set(fig, 'Units', 'Points', 'Position', [100, 100, sizex*nc, sizey*nr])
end
%% Material properties
load(variables, 'mpara');
Emin = mpara(1); Emax = mpara(2); nu = mpara(3);
%% Prepare FEA
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nelem = size(T,1);
nodes = size(X,1);
ndof = 2*nodes;
edofMat = zeros(nelem,8);
edofMat(:,1:2) = [2*T(:,1)-1 2*T(:,1)];
edofMat(:,3:4) = [2*T(:,2)-1 2*T(:,2)];
edofMat(:,5:6) = [2*T(:,3)-1 2*T(:,3)];
edofMat(:,7:8) = [2*T(:,4)-1 2*T(:,4)];
iK = reshape(kron(edofMat,ones(8,1))',64*nelem,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelem,1);
% iF = reshape(edofMat',8*nelem,1);
U = zeros(ndof,1);
% For stress computations
B = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    0 -1/2 0 -1/2 0 1/2 0 1/2
    -1/2 -1/2 -1/2 1/2 1/2 1/2 1/2 -1/2];           % strain-displacement matrix
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       % constitutive matrix - plane stress
dSIGdu = D*B;                                       % for computing adjoint loads
% For stress stifffness matrix
BNL = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    -1/2 0 -1/2 0 1/2 0 1/2 0
    0 -1/2 0 1/2 0 1/2 0 -1/2
    0 -1/2 0 -1/2 0 1/2 0 1/2];          % BNL matrix
%% Prepare augmented PDE filter (with Robin BC)
load(variables, 'rmin', 'xi');
l_o = rmin/2/sqrt(3);             % bulk length scale parameter
l_s = l_o*xi;                     % surface length scale parameter
% PDE filter stiffness matrix
edofMatF = [(edofMat(:,1)+1)/2 (edofMat(:,3)+1)/2 (edofMat(:,5)+1)/2 (edofMat(:,7)+1)/2];
KEF = l_o^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
    helem^2*[4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelem,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelem,1);
sKF = reshape(KEF(:)*ones(1,nelem),16*nelem,1);
% PDE filter "mass" matrix
sMF = 0*sKF;
ME1 = l_s*helem*[2 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 2]/6;                        % mass matrix left   face (1&4)
ME2 = l_s*helem*[2 1 0 0; 1 2 0 0; 0 0 0 0; 0 0 0 0]/6;                        % mass matrix bottom face (1&2)
ME31 = l_s*helem*[0 0 0 0; 0 2 1 0; 0 1 2 0; 0 0 0 0]/6;                       % mass matrix right  face (2&3)
ME41 = l_s*helem*[0 0 0 0; 0 0 0 0; 0 0 2 1; 0 0 1 2]/6;                       % mass matrix top    face (3&4)
% Populate "mass" matrix
idv = zeros(1,nelem);                               % identifies elements for the boundary surface terms
idv = idv*0; idv(T(:,9)==1) = 1;                    % left face
sMF = sMF + reshape(ME1(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,10)==1) = 1;                   % bottom face
sMF = sMF + reshape(ME2(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,11)==1) = 1; idv(solids)=0;    % right face
sMF = sMF + reshape(ME31(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,12)==1) = 1; idv(solids)=0;    % top face
sMF = sMF + reshape(ME41(:)*idv,16*nelem,1);
KF = sparse(iKF,jKF,sKF+sMF);
[LF, FLAG, PF] = chol(KF, 'lower', 'matrix');       % Using permutation matrix 
% Transformation Matrix
iTF = reshape(edofMatF,4*nelem,1);
jTF = reshape(repmat(1:nelem,4,1)',4*nelem,1);
sTF = repmat(1/4,4*nelem,1)*helem^2;
TF = sparse(iTF,jTF,sTF);
filter = @(x) (TF'*(PF*(LF'\(LF\(PF'*(TF*x(:)))))))/helem^2;
clear('edofMatF', 'iKF', 'jKF', 'sKF', 'sMF', ...
    'KF', 'idv', 'iTF', 'jTF', 'sTF')

% BC needed for msolveq and meigen
bc = (1:ndof)';
bc(freedofs) = [];
bc(:, 2) = 0;
%% Initialize MMA
% Initial guess
load(variables, 'x0strat');
switch lower(x0strat)
    case 'homg'
        x = 0.50*ones(nelem, 1);
    case 'comp'
        load(compref, 'x')
    case 'rand'
        x = rand([nelem, 1]);
    otherwise
        x = 0.50*ones(nelem, 1);
end

xold1           = x;                    % xval, one iteration ago (provided that iter>1).
xold2           = x;                    % xval, two iterations ago (provided that iter>2).
m               = 2;                    % number of general constraints.
n               = nelem;                % number of design variables x_j.
mmalen          = 0.200;                % maximum change every design update
xmin            = 0.00*ones(n,1);       % column vector with the lower bounds for the variables x_j.
xmax            = 1e+0*ones(n,1);       % column vector with the upper bounds for the variables x_j.
low             = 0*ones(n,1);          % column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp             = ones(n,1);            % column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0              = 1;                    % the constants a_0 in the term a_0*z.
a               = zeros(m,1);           % column vector with the constants a_i in the terms a_i*z.
c_MMA           = 1000*ones(m,1);       % column vector with the constants c_i in the terms c_i*y_i.
d               = ones(m,1);            % column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
fval            = ones(m,1);            % initial constraint values
chg_phs         = inf;
chg_mma         = inf;
chg_rms         = inf;
loop_chol       = loop;


% Compute initial density field
xTilde = filter(x);
xPhys = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/...
    (tanh(beta*0.5)+tanh(beta*(1-0.5)));
%% Start iteration
while ((loop < minloop ...
        || CONTINUATION_ONGOING ...
        || chg_phs > 1e-2 ...
        ) || fval(1,1) > 1e-3 || fval(2,1) > 1e-3) ...
        && loop < maxloop
    %% Continuation
    loop = loop + 1;
    CONTINUATION_ONGOING = ...
        pE < pphysend...
        || pN < pNend ...
        || beta < betaend;
    CONTINUATION_UPDATE = ...
        CONTINUATION_ONGOING ...
        && loop >= contnext;
    if CONTINUATION_UPDATE
        contnext = loop + pace;
        if pE == pphysend && pN == pNend
            beta = min(betaend, beta + dbeta);
        end
        pN = min(pNend,     pN + dpN);
        pE = min(pphysend,  pE + dpphys);
        pS = min(pphysend,  pS + dpphys);
    end
    %% Conditions for CA solve
    SOLVE_EXACTLY = ...
        loop < CA_START ...
        || ~mod(loop-loop_chol, CHOL_UPDATE_PACE) ...
        || fct_ang < CA_ANG_TOL ...
        || CONTINUATION_UPDATE;
    %% Solve static equation
    r = Emin+xPhys.^pE*(Emax-Emin);
    sK = reshape(KE(:)*r',64*nelem,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    Kf = K(freedofs, freedofs);

    tstart_stat = tic;
    if SOLVE_EXACTLY || SOLVE_STAT_EXACTLY
        loop_chol = loop;
        % SOLVE STATICS WITH CHOLESKY FACTORIZATION
        tstart_fct = tic;
        [R, FLAG, P] = chol(Kf, 'matrix');
        tsol(loop, 1) = toc(tstart_fct);

        U = msolveq(K, F, bc, R, P);
        % Update 'old' quantities
        Rold = R;
        Pold = P;
        xPhysold = xPhys;
    else
        % SOLVE STATICS WITH CA
        % Compute change in stiffness matrix
        dr = (xPhys.^pE - xPhysold.^pE)*(Emax-Emin);
        sdk = reshape(KE(:)*dr',64*nelem,1);
        dK = sparse(iK, jK, sdk); dK = (dK + dK')/2;
        % Solve using CA
        U = msolveq(K, F, bc, Rold, Pold, dK, NBASIS_STAT);
    end
    tsol(loop, 2) = toc(tstart_stat);

    % Compute residuals
    res_stat(loop, 1) = norm(Kf*U(freedofs) - F(freedofs));
    mag_stat(loop, 1) = norm(F(freedofs), 2);
    %% Compliance and its sensitivity
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    comp = sum(sum((Emin+xPhys.^pE*(Emax-Emin)).*ce));
    dc = -pE*(Emax-Emin)*xPhys.^(pE-1).*ce;
    %% Volume and its sensitivity
    v = sum(xPhys)*helem^2/(sizey*sizex);
    dv = ones(size(xPhys))*helem^2/(sizey*sizex);
    %% Compute geometric stiffness
    EPS = U(edofMat)*B';                                % strain
    SIG = EPS*D;                                        % stress for E=1
    sG = 0*sK; dsGdx = sG;
    for el = 1:nelem
        tau = [SIG(el,1) SIG(el,3) 0 0;
            SIG(el,3) SIG(el,2) 0 0 ;
            0 0 SIG(el,1) SIG(el,3);
            0 0 SIG(el,3) SIG(el,2)];
        BNLTtau = BNL'*tau;
        BNLTtauBNL = BNLTtau*BNL;
        l1 = (el-1)*64+1; l2 = el*64;
        sG(l1:l2) = Emax*xPhys(el,1)^pS*helem^2*BNLTtauBNL(:);
        dsGdx(l1:l2) = Emax*pS*xPhys(el,1)^(pS-1)*helem^2*BNLTtauBNL(:);
    end
    KNL = sparse(iK,jK,sG); KNL = (KNL+KNL')/2;
    %% Solve eigenvalue problem
    % COMPUTE EIGS
    tstart_eigs = tic;
    if SOLVE_EXACTLY || SOLVE_EIGS_EXACTLY
        % SOLVE EIGS WITH CHOLESKY FACTORIZATION
        if ~SOLVE_EXACTLY && ~SOLVE_STAT_EXACTLY
            % If static problem was solved using CA
            % we need a fresh factorization of K to solve adjoint exactly
            [R, FLAG, P] = chol(Kf, 'matrix');
        end
        [evecs, evals] = meigenSM(-KNL, K, bc, nevals, R, P);
        PHI = evecs./sqrt(dot(evecs, K*evecs));
        PHIold = PHI;
        % Update CA's orthovecs if 'old' option is used
        CAorthovecs = PHIold(freedofs, :);
    else
        % SOLVE EIGS WITH CA
%         NB = [NBASIS_EIGS*ones(2, 1); 3*ones(nevals-2, 1)];
        NB = NBASIS_EIGS;
        [evecs, evals] = meigenSM(-KNL, K, bc, nevals, ...
            Rold, Pold, dK, PHIold, NB, CAorthotype, CAorthovecs);
        PHI = evecs./sqrt(dot(evecs, K*evecs));
    end
    tsol(loop, 3) = toc(tstart_eigs);

    % Compute residuals
    KNLf = KNL(freedofs, freedofs);
    PHIf = PHI(freedofs, :);
    for k = 1:nevals
        res_eigs(loop, k) = norm(-KNLf*PHIf(:, k) - evals(k, k)*Kf*PHIf(:, k), 2);
        mag_eigs(loop, k) = norm(Kf*PHIf(:, k), 2);
    end

    % Compute blf and pnorm of blf
    mu = diag(evals);
    mu_max_acc = max(mu);
    mu_max_app = norm(mu/max(mu),pN)*max(mu);
    lambda = 1./mu;
    lambda_min_acc = 1/mu_max_acc;
    lambda_min_app = 1/mu_max_app;
    %% Sensitivities of buckling load factor
    % Compute sensitivity term #1 -PHI'*(dGdx+mu*dKdx)*PHI
    dmu = zeros(nelem,nevals);
    for el = 1:nelem
        l1 = (el-1)*64+1; l2 = el*64;
        dsGe = reshape(dsGdx(l1:l2),8,8);
        dsKe = pE*(Emax-Emin)*xPhys(el,1)^(pE-1)*KE;
        dmu11 = sum(PHI(edofMat(el,:)',:).*(dsGe*PHI(edofMat(el,:)',:)));
        dmu12 = mu'.*sum(PHI(edofMat(el,:)',:).*(dsKe*PHI(edofMat(el,:)',:)));
        dmu(el,:) = - (dmu11 + dmu12);
    end
    % Compute adjoint loads
    ADJload = 0*PHI;
    for el = 1:nelem
        edof = edofMat(el,:);
        GradPHI=BNL*PHI(edof,:);
        BNLTdtauBNL = [GradPHI(1,:).^2+GradPHI(3,:).^2;
            GradPHI(2,:).^2+GradPHI(4,:).^2;
            2*(GradPHI(1,:).*GradPHI(2,:)+GradPHI(3,:).*GradPHI(4,:))];
        vals = dSIGdu'*Emax*xPhys(el,1)^pS*helem^2*BNLTdtauBNL;
        ADJload(edof,:) = ADJload(edof,:) - vals;
    end
    ADJsol = 0*ADJload;

    % SOLVE ADJOINTS
    tstart_adjt = tic;
    if SOLVE_EXACTLY || SOLVE_ADJT_EXACTLY
        % SOLVE ADJOINTS WITH CHOLESKY FACTORIZATION
        if ~SOLVE_EXACTLY && ~SOLVE_STAT_EXACTLY && ~SOLVE_EIGS_EXACTLY
            % If static problem and eigenproblem was solved using CA
            % we need a fresh factorization of K to solve adjoint
            [R, FLAG, P] = chol(K(freedofs, freedofs), 'matrix');
        end
        for k = 1:nevals
            ADJsol(:, k) = msolveq(K, ADJload(:, k), bc, R, P);
        end
    else
        % SOLVE ADJOINTS WITH CA
        for k = 1:nevals
%             NBK = NBASIS_ADJT*(k <= 2) + 3*(k > 2);
            NBK = NBASIS_ADJT;
            ADJsol(:, k) = msolveq(K, ADJload(:, k), bc, ...
                Rold, Pold, dK, NBK);
        end
    end
    tsol(loop, 4) = toc(tstart_adjt);

    % Compute sensitivities
    for j = 1:nevals
        adjsol = ADJsol(:,j);
        vals = sum((adjsol(edofMat)*KE).*U(edofMat),2);
        dmu(:,j) = dmu(:,j)-pE*(Emax-Emin)*xPhys.^(pE-1).*vals;
    end
    dmucdmu = (mu/mu_max_app)'.^(pN-1);
    dmu_max_app = dmu*dmucdmu';
    
    % Compute residuals
    ADJsolf = ADJsol(freedofs, :);
    ADJloadf = ADJload(freedofs, :);
    for k = 1:nevals
        res_adjt(loop, k) = norm(Kf*ADJsolf(:, k) - ADJloadf(:, k), 2);
        mag_adjt(loop, k) = norm(ADJloadf(:, k), 2);
    end
    %% Chain rule for projection and filter
    dxPhys = (1 - (tanh(beta*(xTilde-0.5))).^2)*beta/(tanh(beta*0.5)+tanh(beta*(1-0.5)));
    dv          = filter(dv.*dxPhys);
    dc          = filter(dc.*dxPhys);
    dmu_max_app = filter(dmu_max_app.*dxPhys);
    %% Draw design and stress
    if SHOW_DESIGN
        for k = 1:nc*nr
            ck = floor((k-1)/nr)+1;
            rk = mod(k-1, nr)+1;
            % Plot stuff
            if ck == 1
                % First column is xPhys and objective and constraints
                switch rk
                    case 1
                        xx = xPhys;
                        tit = '\nu';
                    case 2
                        xx = log(abs(dmu_max_app));
                        tit = '\partial \mu_c';
                    case 3
                        xx = dc;
                        tit = '\partial c';
                end
            else
                xx = log(abs(dmu(:, rk)));
                tit = sprintf('\\partial \\mu_%i', rk);
            end
            v_img = 2*(xx-min(xx))/(max(xx)-min(xx)) - 1;
            imgk = sparse(i_img,j_img,v_img);
            axes(ax(k));                                        %#ok<LAXES>
            imagesc(imgk);

            % Title
            xt = ruler2num(sizex/helem*0.50, ax(k).XAxis);
            yt = ruler2num(sizey/helem*1.00, ax(k).YAxis);
            text(ax(k), xt, yt, tit, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top');

            pause(0.05);                    % To make sure plots are right size and in order
        end
        drawnow;
    end
    %% MMA
    xval  = x;
    if (loop==1)
        scale = 10/comp;
    end
    f0val = scale*comp;
    df0dx = scale*dc;
    fval(1,1) = v/volfrac - 1;
    dfdx(1,:) = dv/volfrac;
    fval(2,1) = mu_max_app/mustar - 1;
    dfdx(2,:) = dmu_max_app/mustar;

    xmink = max(xmin,xval-mmalen);
    xmaxk = min(xmax,xval+mmalen);
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,xval,xmink,xmaxk,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    %% Update MMA Variables
    x        = xmma;
    xold2    = xold1;
    xold1    = xval;
    xPhysold1 = xPhys;

    % Filter new fields
    xTilde      = filter(x);
    xPhys       = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/(tanh(beta*0.5)+tanh(beta*(1-0.5)));

    % Compute design changes
    dx = abs(x - xold1);
    dxPhys = abs(xPhys - xPhysold1);

    chg_mma     = max(dx);
    chg_phs     = max(dxPhys);
    chg_rms     = rms(dxPhys);
    chg_ang     = dot(xPhys/norm(xPhys, 2), xPhysold1/norm(xPhysold1, 2));
    fct_ang     = dot(xPhys/norm(xPhys, 2), xPhysold/norm(xPhysold, 2));
    %% Print results
    SV = bitor([SOLVE_STAT_EXACTLY, SOLVE_EIGS_EXACTLY, SOLVE_ADJT_EXACTLY], SOLVE_EXACTLY);
    fprintf([...
        'ITER: %3i OBJ: %+10.3e CONST: ', repmat('%+10.3e ', 1, m), ...
        'BLF: ', repmat('%+10.3e ', 1, nevals), 'BLFAPP %+10.3e ', ...
        'CH: %5.3f CHPHS: %5.3f RMSPHS: %5.3f CHGANG: %5.3f FCTANG: %5.3f ', ...
        'pE: %5.2f pN: %5.2f BETAHS: %3i ', ...
        'EXACT(S/E/A): %1i %1i %1i \n'],...
        loop,f0val,fval',lambda',lambda_min_app,chg_mma,chg_phs,chg_rms,chg_ang,fct_ang,pE,pN,beta,SV);
    %% Save data
    stats(loop,1:11+m+nevals) = [f0val,fval',lambda',chg_mma,chg_phs,chg_rms,chg_ang,pE,pN,beta,SV];
end
%% Save data etc
runtime = toc(tstart);
profile_data = profile('info');
% Save figure
if SHOW_DESIGN
    pause(2.00); saveas(gcf, 'design.png'); pause(2.00);
end
% Save results
DOMVARS = load(domaindata);
INPVARS = load(variables);
save('results.mat', ...
    'xPhys', 'xTilde', 'x', 'stats', 'tsol', 'profile_data',...
    'INPVARS', 'DOMVARS', 'T', ...
    'res_stat', 'mag_stat', ...
    'res_eigs', 'mag_eigs', ...
    'res_adjt', 'mag_adjt');

fprintf('\nFinished job at %s \n', datestr(datetime('now')));
