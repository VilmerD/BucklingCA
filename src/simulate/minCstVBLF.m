%% Minimum max mu (1/lambda) s.t. volume
% Approximation of max mu using p-norm
% Single density field with eta=0.5
profile off
profile on -history;
addpath(genpath('~/Documents/Projects/BucklingCA/data/'))
addpath(genpath('~/Documents/Projects/BucklingCA/src/generate_geometry'))
addpath(genpath('~/Documents/Projects/BucklingCA/src/solvers'))
tstart = tic;
%% Load run data
% If there is an inputfile named input.mat then load data from it, otherwise
% use default parameters
if isfile(fullfile(cd, 'input.mat'))
    variables = 'input.mat';
else
    variables = 'data/default_data.mat';
end
%% Optimization parameters
load(variables, 'nevals');                         % number of eigenvalues to consider
savefile = 'results.mat';
%% Domain size and discretization
% element size (all elements are square)
load(variables, 'helem', 'domain');
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
%% Target blf
comp_ref = sprintf('%s.mat', domain);
load(fullfile('data/compliance_reference', comp_ref), 'lambda');
lamstar = lambda(1);        % Critical blf
load(variables, 'sf');      % Target blf is safety factor * critical blf
lamstar = lamstar*sf;
mustar = 1/lamstar;
%% Initialize optimization
minloop = 350;                                      % minimum number of loops
maxloop = 400;

pE = 2;                                             % SIMP penalty for linear stiffness
pS = 2;                                             % SIMP penalty for stress stiffness
pphysmax = 6;                                       % maximum penalization
dpphys = 0.25;                                      % change in penalty
psteps = ceil((pphysmax-pE)/dpphys);

pN = 8;                                             % p-norm for eigenvalues
pNmax = 64;
dpN = (pNmax - pN)/psteps;

beta = 6;                                           % thresholding steepness
betamax = 6;
dbeta = 0;

jumpnext = 30;                                      % next jump loop
pace = (minloop-jumpnext)/psteps;                   % update pace

loop = 0;
%% Allocate some data to track solution
stats = zeros(minloop,10+3+nevals);
tsol = zeros(maxloop, 4);
res_stat = zeros(minloop, 1);
res_eigs = zeros(minloop, nevals);
res_adjt = zeros(minloop, nevals);
%% Initialize CA
% Booleans controlling whether or not to use CA for
% the individual problems
SOLVE_STAT_EXACTLY = 0;
SOLVE_EIGS_EXACTLY = 0;
SOLVE_ADJT_EXACTLY = 0;
% Number of basis vectors
load(variables, 'NBASIS_STAT', 'NBASIS_EIGS', 'NBASIS_ADJT');
% Orthogonalization type for eigenproblem
CAOPTS.orthotype = 'old';
CAOPTS.orthovecs = [];          % Old eigenmodes if orthotype is 'old', else empty
% When CA should be started, and how often to factorize
load(variables, 'CA_START', 'CHOL_UPDATE_PACE', 'CA_CHANGE_TOL')
%% Initialize figure
nc = 2;
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
                'Position', [(c-1)/nc, (nr-r)/nr, 1/nc, 1/nr]);      %#ok<SAGROW>
            hold(ax(i), 'on');
            axis(ax(i), 'off');
            axis(ax(i), 'tight');
            colormap(ax(i), 'winter');
            pause(0.05);
        end
    end
    set(fig, 'Units', 'Points', 'Position', [100, 100, sizex*nc, sizey*nr])
end
%% Material properties
Emax = 2e5;
Emin = Emax*1e-6;
nu = 0.3;
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
load(variables, 'rmin');
xi = 1;                     % default for Robin BC (ratio l_s/l_o)
xi_corner = xi;             	% large xi near re-entrant corner
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
[LF, FLAG, PF] = chol(KF,'lower', 'matrix');
% Transformation Matrix
iTF = reshape(edofMatF,4*nelem,1);
jTF = reshape(repmat(1:nelem,4,1)',4*nelem,1);
sTF = repmat(1/4,4*nelem,1)*helem^2;
TF = sparse(iTF,jTF,sTF);
filter = @(x) (TF'*(PF*(LF'\(LF\(PF'*(TF*x(:)))))))/helem^2;            % actual integration of NdA

% BC needed for msolveq and meigen
bc = (1:ndof)';
bc(freedofs) = [];
bc(:, 2) = 0;
%% Initialize MMA
m               = 2;                    % number of general constraints.
n               = nelem;                % number of design variables x_j.
mmalen          = 0.100;                % maximum change every design update
xmin            = 1e-6*ones(n,1);       % column vector with the lower bounds for the variables x_j.
xmin(solids)    = 1-1e-3;               % lower bound for solids
xmax            = 1e+0*ones(n,1);       % olumn vector with the upper bounds for the variables x_j.
xmax(voids)     = 1e-3;                 % upper bound for voids
x               = 0.55*ones(n, 1);
xold1           = x;                    % xval, one iteration ago (provided that iter>1).
xold2           = x;                    % xval, two iterations ago (provided that iter>2).
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
loop_chol       = 1;

% Compute initial density field
xTilde = filter(x);
xPhys = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/...
    (tanh(beta*0.5)+tanh(beta*(1-0.5)));
%% Start iteration
while (((loop < jumpnext || CONTINUATION_ONGOING) ...
        || pE < pphysmax ...
        || chg_phs > 5e-2 ...
        ) || fval(1,1) > 1e-3) && loop < maxloop
    %% Continuation
    loop = loop + 1;
    CONTINUATION_ONGOING = ...
        pE < pphysmax...
        || pN < pNmax ...
        || beta < betamax;
    CONTINUATION_UPDATE = ...
        CONTINUATION_ONGOING ...
        && loop >= jumpnext;
    if CONTINUATION_UPDATE
        jumpnext = loop + pace;
        if pE == pphysmax
            beta = min(betamax, beta + dbeta);
        end
        pN = min(pNmax,     pN + dpN);
        pE = min(pphysmax,  pE + dpphys);
        pS = min(pphysmax,  pS + dpphys);
    end
    %% Conditions for CA solve
    SOLVE_EXACTLY = ...
        loop < CA_START ...
        || ~mod(loop-loop_chol, CHOL_UPDATE_PACE) ...
        || chg_phs > CA_CHANGE_TOL ...
        || CONTINUATION_UPDATE;
    %% Solve static equation
    sK = reshape(KE(:)*(Emin+(xPhys').^pE*(Emax-Emin)),64*nelem,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    tstart_stat = tic;
    if SOLVE_EXACTLY || SOLVE_STAT_EXACTLY
        loop_chol = loop;
        % SOLVE STATICS WITH CHOLESKY FACTORIZATION
        tstart_fct = tic;
        [R, FLAG, P] = chol(K(freedofs, freedofs), 'matrix');
        tsol(loop, 1) = toc(tstart_fct);
        U = msolveq(K, F, bc, R, P);
        % Update 'old' quantities
        Rold = R;
        Pold = P;
        Kold = K;
    else
        % SOLVE STATICS WITH CA
        U = msolveq(K, F, bc, ...
            Rold, Pold, Kold, NBASIS_STAT);
    end
    tsol(loop, 2) = toc(tstart_stat) - tsol(loop, 1);

    % Compute residuals
    res_stat(loop, 1) = norm(K(freedofs, freedofs)*U(freedofs) - F(freedofs));
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
        NB = [NBASIS_EIGS*ones(2, 1); 2*ones(nevals-2, 1)];
        [evecs, evals] = meigenSM(-KNL, K, bc, nevals, ...
            Rold, Pold, Kold, PHIold, NB, CAOPTS);
        PHI = evecs./sqrt(dot(evecs, K*evecs));
    end
    tsol(loop, 3) = toc(tstart_eigs);

    % Compute residuals
    for k = 1:nevals
        res_eigs(loop, k) = norm(KNL(freedofs, freedofs)*PHI(freedofs, k) - evals(k)*K(freedofs, freedofs)*PHI(freedofs, k));
        res_eigs(loop, k) = res_eigs(loop, k)/norm(K(freedofs, freedofs)*PHI(freedofs, k));
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
            NBK = NBASIS_ADJT*(k <= 2) + 2*(k > 2);
            ADJsol(:, k) = msolveq(K, ADJload(:, k), bc, ...
                Rold, Pold, Kold, NBK);
        end
    end
    tsol(loop, 4) = toc(tstart_adjt);
    % Compute residuals
    for k = 1:nevals
        res_adjt(loop, k) = norm(K*ADJsol(:, k) - ADJload(:, k));
    end

    % Compute sensitivities
    dmu2 = 0*dmu1;
    for j = 1:nevals
        adjsol = ADJsol(:,j);
        vals = sum((adjsol(edofMat)*KE).*U(edofMat),2);
        dmu2(:,j) = -pE*(Emax-Emin)*xPhys.^(pE-1).*vals;
    end
    dmu = dmu1 + dmu2;
    dmucdmu = (mu/mu_max_app)'.^(pN-1);
    dmu_max_app = dmu*dmucdmu';
    %% Chain rule for projection and filter
    dxPhys = (1 - (tanh(beta*(xTilde-0.5))).^2)*beta/(tanh(beta*0.5)+tanh(beta*(1-0.5)));
    dv          = filter(dv.*dxPhys);
    dc          = filter(dc.*dxPhys);
    dmu_max_app = filter(dmu_max_app.*dxPhys);
    %% Draw design and stress
    if SHOW_DESIGN
        for k = 1:nc*nr
            % Plot stuff
            switch k
                case 1
                    xx = xPhys;
                    tit = '\nu';
                case 2
                    xx = log(abs(dmu_max_app));
                    tit = '\partial \mu_c';
                case 3
                    xx = dc;
                    tit = '\partial c';
                otherwise
                    xx = log(abs(dmu(:, k-nr)));
                    tit = sprintf('\\partial \\mu_%i', k-nr);
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
            pause(0.25);    % To make sure plots are right size and in order
        end
        drawnow;
    end
    %% MMA
    xval  = x;
    if (loop==1)
        scale = 100/comp;
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
    xnew     = xmma;
    xold2    = xold1;
    xold1    = xval;

    % Filter new fields
    xTildenew   = filter(xnew);
    xPhysnew    = (tanh(beta*0.5)+tanh(beta*(xTildenew-0.5)))/(tanh(beta*0.5)+tanh(beta*(1-0.5)));
    chg_mma     = max(abs(xnew-xval));
    chg_phs     = max(abs(xPhysnew - xPhys));
    chg_rms     = rms(xPhysnew - xPhys);

    % Update fields
    x = xnew;
    xTilde = xTildenew;
    xPhys = xPhysnew;
    %% Print results
    SV = bitor([SOLVE_STAT_EXACTLY, SOLVE_EIGS_EXACTLY, SOLVE_ADJT_EXACTLY], SOLVE_EXACTLY);
    fprintf([...
        'ITER: %3i OBJ: %+10.3e CONST: ', repmat('%+10.3e ', 1, m), ...
        'BLF: ', repmat('%+10.3e ', 1, nevals), 'BLFAPP %+10.3e ', ...
        'CH: %5.3f CHPHS: %5.3f RMSPHS: %5.3f ', ...
        'pE: %5.2f pN: %5.2f BETAHS: %3i ', ...
        'EXACT(S/E/A): %1i %1i %1i \n'],...
        loop,f0val,fval',lambda',lambda_min_app,chg_mma,chg_phs,chg_rms,pE,pN,beta,SV);
    %% Save data
    stats(loop,1:10+m+nevals) = [f0val,fval',lambda',chg_mma,chg_phs,chg_rms,pE,pN,beta,SV];
end
runtime = toc(tstart);
profile_data = profile('info');
saveas(gcf, 'design.png');
close(fig)
clearvars i* j* s* d* edof* *old* *new* *mma* PHI* ...
    KF TF LF PF filter ...
    K R P KNL ...
    evecs adj* ADJ* EPS SIG ...
    ax low upp ...
    freedofs X U ce *img* *val* ...
    -EXCEPT xPhys xTilde x sizex sizey domain T stats
save('results.mat');
