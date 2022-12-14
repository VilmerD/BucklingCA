function [X, L, varargout] = meigenSM(A, B, bc, ne, Rold, Pold, Bold, Xold, s, options)
% [X, L] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem (A+l*B)x = 0 with homogeneous boundary 
% conditions bc
%
% [X, L, R, P] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem and returns R, the cholesky factorization 
% of the free part of K
%
% [X, L, d, B] = MEIGENSM(A, B, bc, ne, Rold, Bold, Xold, s, options)
% finds the first ne eigenvalues using CA. Kold is the stiffness matrix
% corresponding to the factorization R and the eigenmodes psi0. CA uses s
% basis vectors to compute the eigenvalues, and options is a struct
% containing information for the CA procedure. The basis vectors for the
% reduced order model are computed using CAeigs.
%

% Extracting free part of vibrational problem
dofsfree = (1:size(A, 1))';
dofsfree(bc(:, 1)) = [];

Aff = A(dofsfree, dofsfree);
Bff = B(dofsfree, dofsfree);

% Solving generalized eigenvalue problem
if nargin <= 6
    % Solve the problem using cholesky factorization
    if nargin < 6
        [Rcurr, ~, Pcurr] = chol(Aff, 'matrix');
    else
        Rcurr = Rold;
        Pcurr = Pold;
    end
    [pcurr, ~, ~] = find(Pcurr);
    [Pf, L] = eigs(Aff, Rcurr, ne, 'largestreal', ...
        'IsCholesky', true, 'CholeskyPermutation', pcurr, ...
        'Display', false);
    
    % Normalize with respect to the mass matrix
    Pf = Pf./sqrt(dot(Pf, Bff*Pf));
    
    % Typically the user wants the cholesky factorization
    varargout = {Rcurr, Pcurr};
else
    % Solve the problem using a reduced order model if the proper
    [Pf, L, d, B] = CAeigs(Aff, Bff, ne, Rold, Pold, ...
        Bold(dofsfree, dofsfree), Xold(dofsfree, :), s, options);
    varargout = {d, B};
end

% Sorting eigenvalues in ascending order
L = reshape(diag(L), ne, 1);
[Ld, I] = sort(L, 'descend');
L = diag(Ld);
Pf = Pf(:, I);

% Making eigenvectors full
X = zeros(size(A, 1), ne);
X(dofsfree, :) = Pf;
end