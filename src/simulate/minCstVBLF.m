%% Minimum max gamma (1/lambda) s.t. volume
% Approximation of max gamma using p-norm
% Single density field with eta=0.5
clear;
% close all;
fig = figure(1);
commandwindow;
clc;
tic;
%% Choose problem type
% prtype = 'WEIGHTED';   % Gradually adds weight to x'*(1-x) to promote discrete design
% prtype = 'BLF';                   % Maximize the buckling load factor
% prtype =  'COMP';                 % Minimize compliance
prtype =  'COMP_ST_BLF';            % Minimize compliance subject to buckling const.
weight = 1;
mustar = 1/10;
%% Domain size and discretization
domain = 'twobar';                  % options are: 'column' 'spire' 'twobar'
sizex = 150;                         % physical size in x-direction
sizey = 300;                        % physical size in y-direction
helem = 001;                          % element size (all elements are square)
fig.Position(3:4) = [sizex*1.20 sizey*3*1.20]/sizex*400;
%% Optimization parameters
pE = 3;                             % SIMP penalty for linear stiffness
pS = 3;                             % SIMP penalty for stress stiffness
pN = 8;                             % p-norm for eigenvalues
rmin = 5;                           % filter radius (for convolution density filter)
nevals = 6;                         % number of eigenvalues to consider
filename = sprintf('processed_data/%s_trial%i.mat', domain, numel(dir('processed_data'))-1);
%% Material properties
Emax = 2e5;
Emin = Emax*1e-6;
nu = 0.3;
%% Prepare design domain
% X = nodal coordinates
% T = element connectivity and data
% i_img,j_img = for displaying design as a matrix using imagesc
switch domain
    case 'column'
        [X,T,i_img,j_img,solids,voids,F,freedofs] = generate_column(sizex,sizey,helem,false);
    case 'spire'
        [X,T,i_img,j_img,solids,voids,F,freedofs] = generate_spire(sizex,sizey,helem,false);
    case 'twobar'
        [X,T,i_img,j_img,solids,voids,F,freedofs] = generate_twobar(sizex,sizey,helem,false);
end
%% Prepare FEA (88-line style)
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
iF = reshape(edofMat',8*nelem,1);
U = zeros(ndof,1); PHI = zeros(ndof,nevals);
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
xi = 1;                     % default for Robin BC (ratio l_s/l_o)
xi_corner = 1;             	% large xi near re-entrant corner
l_o = rmin/2/sqrt(3);       % bulk length scale parameter
l_s = l_o*xi;               % default surface length scale parameter
% PDE filter stiffness matrix
edofMatF = [(edofMat(:,1)+1)/2 (edofMat(:,3)+1)/2 (edofMat(:,5)+1)/2 (edofMat(:,7)+1)/2];
KEF = l_o^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
    helem^2*[4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelem,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelem,1);
sKF = reshape(KEF(:)*ones(1,nelem),16*nelem,1);
% PDE filter "mass" matrix
sMF = 0*sKF;
ME1 = l_s*helem*[2 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 2]/6;                        % mass matrix left face (1&4)
ME2 = l_s*helem*[2 1 0 0; 1 2 0 0; 0 0 0 0; 0 0 0 0]/6;                        % mass matrix bottom face (1&2)
ME31 = l_s*helem*[0 0 0 0; 0 2 1 0; 0 1 2 0; 0 0 0 0]/6;                       % mass matrix right face, standard (2&3)
ME32 = l_s*helem*(xi_corner/xi)*[0 0 0 0; 0 2 1 0; 0 1 2 0; 0 0 0 0]/6;        % mass matrix right face, corner (2&3)
ME41 = l_s*helem*[0 0 0 0; 0 0 0 0; 0 0 2 1; 0 0 1 2]/6;                       % mass matrix top face, standard (3&4)
ME42 = l_s*helem*(xi_corner/xi)*[0 0 0 0; 0 0 0 0; 0 0 2 1; 0 0 1 2]/6;        % mass matrix top face, corner (3&4)
% Populate "mass" matrix
idv = zeros(1,nelem);                               % identifies elements for the boundary surface terms
idv = idv*0; idv(T(:,9)==1) = 1;                    % left face
sMF = sMF + reshape(ME1(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,10)==1) = 1;                   % bottom face
sMF = sMF + reshape(ME2(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,11)==1) = 1; idv(solids)=0;    % right face standard
sMF = sMF + reshape(ME31(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,11)==2) = 1; idv(solids)=0;    % right face corner
sMF = sMF + reshape(ME32(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,12)==1) = 1; idv(solids)=0;    % top face standard
sMF = sMF + reshape(ME41(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,12)==2) = 1; idv(solids)=0;    % top face corner
sMF = sMF + reshape(ME42(:)*idv,16*nelem,1);
KF = sparse(iKF,jKF,sKF+sMF);
LF = chol(KF,'lower');
% Transformation Matrix
iTF = reshape(edofMatF,4*nelem,1);
jTF = reshape(repmat(1:nelem,4,1)',4*nelem,1);
sTF = repmat(1/4,4*nelem,1)*helem^2;
TF = sparse(iTF,jTF,sTF);                           % actual integration of NdA
TF0 = sparse(iTF, jTF, sTF/helem^2); 				% used to average nodal values
%% Initialize optimization
maxloop = 300;
volfrac = 0.16;
x = volfrac*ones(nelem,1);
x(T(:,5)==1) = volfrac*(1e-6);    % voids
x(T(:,5)==2) = volfrac*(1-1e-6);  % solids
beta = 1;
njumps = 8;
betamax = 16;
dbeta = (betamax/beta)^(1/njumps); % factor for multiplying beta
dweight = 1;
weightmin = 1;
changetol = 5e-2;
jumpnext = 60;
pace = min(20,maxloop/(njumps+1));
loop = 0;
Stats = zeros(maxloop,20);
%% Initialize CA
% Booleans controlling whether or not to use CA for
% the individual problems
SOLVE_STAT_EXACTLY = 1;
SOLVE_EIGS_EXACTLY = 1;
SOLVE_ADJT_EXACTLY = 1;
% Number of basis vectors
NBASIS_STAT = 08;
NBASIS_EIGS = 04;
NBASIS_ADJT = 10;
% Orthogonalization type for eigenproblem
CAOPTS.orthotype = 'current';
CAOPTS.orthovecs = [];          % Old eigenmodes if orthotype is 'old', else empty
% When CA should be started, and how often to factorize
CA_START = 20;
CHOL_UPDATE_PACE = 5;

% BC needed for msolveq and meigen
bc = (1:ndof)';
bc(freedofs) = [];
bc(:, 2) = 0;
%% Initialize MMA
m     = 1 + 1*strcmp(prtype,'COMP_ST_BLF');     % number of general constraints.
n     = nelem;                                  % number of design variables x_j.
xmin  = 1e-6*ones(n,1);     % column vector with the lower bounds for the variables x_j.
xmin(solids) = 1-1e-3;      % lower bound for solids
xmax  = ones(n,1);          % olumn vector with the upper bounds for the variables x_j.
xmax(voids) = 1e-3;         % upper bound for voids
xold1 = x;                  % xval, one iteration ago (provided that iter>1).
xold2 = x;                  % xval, two iterations ago (provided that iter>2).
low   = 0*ones(n,1);        % column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);          % column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                  % the constants a_0 in the term a_0*z.
a     = zeros(m,1);         % column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);     % column vector with the constants c_i in the terms c_i*y_i.
d     = ones(m,1);          % column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
fval = ones(m,1);           % initial constraint values

% Compute density field
xTilde = (TF0'*(LF'\(LF\(TF*x(:)))));
xPhys = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/...
    (tanh(beta*0.5)+tanh(beta*(1-0.5)));
%% Start iteration
while ((loop < maxloop || beta < betamax || weight > weightmin) || fval(1,1) > 1e-3)
    %% Continuation
    CONTINUATION_UPDATE = ...
        loop >= jumpnext && ...
        change_phys <= changetol;
    if CONTINUATION_UPDATE
        if beta == betamax
            weight = max(weight*dweight,weightmin);
        end
        beta = min(beta*dbeta,betamax);
        jumpnext = jumpnext + pace;
    end
    loop = loop + 1;
    %% Conditions for CA solve
    SOLVE_EXACTLY = ...
        loop < CA_START ...
        || ~mod(loop, CHOL_UPDATE_PACE) ...
        || CONTINUATION_UPDATE;
    %% Solve static equation
    sK = reshape(KE(:)*(Emin+(xPhys').^pE*(Emax-Emin)),64*nelem,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    if SOLVE_EXACTLY || SOLVE_STAT_EXACTLY
        % SOLVE STATICS WITH CHOLESKY FACTORIZATION
        [R, FLAG, P] = chol(K(freedofs, freedofs), 'matrix');
        [U, ~, ~] = msolveq(K, F, bc, R, P);
        % Update 'old' quantities
        Rold = R;
        Pold = P;
        Kold = K;
    else
        % SOLVE STATICS WITH CA
        U = msolveq(K, F, bc, ...
            Rold, Pold, Kold, NBASIS_STAT);
    end
    %% Compliance and its sensitivity
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    comp = sum(sum((Emin+xPhys.^pE*(Emax-Emin)).*ce));
    dc = -pE*(Emax-Emin)*xPhys.^(pE-1).*ce;
    %% Volume and its sensitivity
    v = sum(xPhys)*helem^2/(sizey*sizex);
    dv = ones(size(xPhys))*helem^2/(sizey*sizex);
    %% MND and its sensitivity
    mnd = (xPhys'*(1-xPhys))*helem^2;
    dmnd = (ones(nelem,1) - 2*xPhys)*helem^2;
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
    %% Solve eigenvalue problem
    KNL = sparse(iK,jK,sG); KNL = (KNL+KNL')/2;
    % COMPUTE EIGS
    if SOLVE_EXACTLY || SOLVE_EIGS_EXACTLY
        % SOLVE EIGS WITH CHOLESKY FACTORIZATION
        if ~SOLVE_EXACTLY && ~SOLVE_STAT_EXACTLY
            % If static problem was solved using CA
            % we need a fresh factorization of K to solve adjoint exactly
            [R, FLAG, P] = chol(K(freedofs, freedofs), 'matrix');
        end
        [evecs, evals] = meigenSM(-KNL, K, bc, nevals, R, P);
        PHI = evecs./sqrt(dot(evecs, K*evecs));
        PHIold = PHI;
        % Update CA's orthovecs if 'old' option is used
        if strcmpi(CAOPTS.orthotype, 'old')
            CAOPTS.orthovecs = PHIold(freedofs, :);
        end
    else
        % SOLVE EIGS WITH CA
        [evecs, evals, eigdelt] = meigenSM(-KNL, K, bc, nevals, ...
            Rold, Pold, Kold, PHIold, NBASIS_EIGS, CAOPTS);
        PHI = evecs./sqrt(dot(evecs, K*evecs));
    end

    mu = diag(evals);
    mu_max_acc = max(mu);
    mu_max_app = norm(mu,pN);
    lambda = 1./mu;
    lambda_min_acc = 1/mu_max_acc;
    lambda_min_app = 1/mu_max_app;
    %% Sensitivities of buckling load factor
    % Compute sensitivity term #1 -PHI'*(dGdx+mu*dKdx)*PHI
    dmu1 = zeros(nelem,nevals);
    for el = 1:nelem
        l1 = (el-1)*64+1; l2 = el*64;
        dsGe = reshape(dsGdx(l1:l2),8,8);
        dsKe = pE*(Emax-Emin)*xPhys(el,1)^(pE-1)*KE;
        dmu11 = sum(PHI(edofMat(el,:)',:).*(dsGe*PHI(edofMat(el,:)',:)));
        dmu12 = mu'.*sum(PHI(edofMat(el,:)',:).*(dsKe*PHI(edofMat(el,:)',:)));
        dmu1(el,:) = - (dmu11 + dmu12);
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
            ADJsol(:, k) = msolveq(K, ADJload(:, k), bc, ...
                Rold, Pold, Kold, NBASIS_ADJT);
        end
    end
    
    dmu2 = 0*dmu1;
    for j = 1:nevals
        adjsol = ADJsol(:,j);
        vals = sum((adjsol(edofMat)*KE).*U(edofMat),2);
        dmu2(:,j) = -pE*(Emax-Emin)*xPhys.^(pE-1).*vals;
    end
    dmu = dmu1 + dmu2;
    term1 = 1/pN*(sum(mu.^pN))^(1/pN-1);
    term2 = (mu.^(pN-1))'*dmu';
    dmu_max_app = term1*pN*term2';
    %% Chain rule for projection and filter
    dxPhys = (1 - (tanh(beta*(xTilde-0.5))).^2)*beta / ...
        (tanh(beta*0.5)+tanh(beta*(1-0.5)));
    dv =                            dv.*dxPhys;
    dc =                            dc.*dxPhys;
    dmu_max_app =                   dmu_max_app.*dxPhys;
    dmnd =                          dmnd.*dxPhys;
    dv =                            TF0'*(LF'\(LF\(TF*dv)));
    dc =                            TF0'*(LF'\(LF\(TF*dc)));
    dmu_max_app =                   TF0'*(LF'\(LF\(TF*dmu_max_app)));
    dmnd =                          TF0'*(LF'\(LF\(TF*dmnd)));
    %% Draw design and stress
%     figure(1);
    clf;
    v_img = xPhys;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,1);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('xPhys');
    v_img = dmu_max_app;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,2);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('d{\mu}_{max}');
    v_img = dc;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,3);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('d{c}');
    drawnow;
    %% MMA
    xval  = x;
    switch prtype
        case 'COMP' % minCstV
            if (loop==1)
                scale = 10/comp;
            end
            f0val = scale*comp;
            df0dx = scale*dc;
            fval(1,1) = v/volfrac - 1;
            dfdx(1,:) = dv/volfrac;
        case 'BLF' % min gamma
            if (loop==1)
                scale = 10/mu_max_app;
            end
            f0val = scale*mu_max_app;
            df0dx = scale*dmu_max_app;
            fval(1,1) = v/volfrac - 1;
            dfdx(1,:) = dv/volfrac;
        case 'WEIGHTED'
            if (loop==1)
                scale1 = 10/mu_max_app;
                scale2 = 10/mnd;
            end
            f0val = scale1*weight*mu_max_app + scale2*(1-weight)*mnd;
            df0dx = scale1*weight*dmu_max_app + scale2*(1-weight)*dmnd;
            fval(1,1) = v/volfrac - 1;
            dfdx(1,:) = dv/volfrac;
        case 'COMP_ST_BLF'
            if (loop==1)
                scale1 = 100/comp;
                scale2 = 10/mnd;
            end
            f0val = scale1*weight*comp + scale2*(1-weight)*mnd;
            df0dx = scale1*weight*dc + scale2*(1-weight)*dmnd;
            fval(1,1) = v/volfrac - 1;
            dfdx(1,:) = dv/volfrac;
            fval(2,1) = mu_max_app/mustar - 1;
            dfdx(2,:) = dmu_max_app/mustar;
    end
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,xval,max(xmin,xval-0.2),min(xmax,xval+0.2),xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % Update MMA Variables
    xnew     = xmma;
    xold2    = xold1;
    xold1    = xval;

    xTildenew = (TF0'*(LF'\(LF\(TF*xnew))));
    xPhysnew = (tanh(beta*0.5)+tanh(beta*(xTildenew-0.5)))/...
        (tanh(beta*0.5)+tanh(beta*(1-0.5)));
    change      = max(abs(xnew-xval));
    change_phys = max(abs(xPhysnew - xPhys));


    x = xnew;
    xTilde = xTildenew;
    xPhys = xPhysnew;

    %% Print results
    SOLVE_EXACTLY_VEC = [...
        SOLVE_STAT_EXACTLY || SOLVE_EXACTLY...
        SOLVE_EIGS_EXACTLY || SOLVE_EXACTLY...
        SOLVE_ADJT_EXACTLY || SOLVE_EXACTLY];
    fprintf([...
        'ITER: %3i OBJ: %+10.3e CONST: ', repmat('%+10.3e ', 1, m), ...
        'CH: %5.3f CHPHS: %5.3f BETAHS: %6.3f WEIGHT: %4.3f', ...
        'EXACT(S/E/A): %1i %1i %1i \n'],...
        loop,f0val,fval',change,change_phys,beta,weight,SOLVE_EXACTLY_VEC(:));
    %% Save data
    Stats(loop,1:1+m+1+nevals+1) = [f0val fval' beta lambda' weight];
end
runtime = toc;
save(filename);



