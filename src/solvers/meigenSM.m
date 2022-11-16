function [X, L, varargout] = meigenSM(A, B, bc, ne, Rold, Pold, dB, Xold, s, ...
    orthotype, orthovecs)
% [X, L] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem (A+l*B)x = 0 with homogeneous boundary 
% conditions bc
%
% [X, L, R, P] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem and returns R, the cholesky factorization 
% of the free part of A
%
% [X, L, d, B] = MEIGENSM(A, B, bc, ne, Rold, dB, Xold, s, options)
% finds the first ne eigenvalues using CA. dB is the change in B
% corresponding to the factorization Rold and the eigenmodes Xold. CA uses s
% basis vectors to compute the eigenvalues. The basis vectors for the
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
    elseif nargin == 5
        % No permutation matrix is given, assume no permutation
        Rcurr = Rold;
        Pcurr = speye(size(Rcurr));
    else
        % Both cholesky and P is given, just rename
        Rcurr = Rold;
        Pcurr = Pold;
    end
    % Since eigs requires 'vector' style of P, but 'matrix' is used
    % elsewhere we must transform, which is easily done using find
    [pcurr, ~, ~] = find(Pcurr);
    % Solve using eigs
    [Pf, L] = eigs(Aff, Rcurr, ne, 'largestreal', ...
        'IsCholesky', true, ....
        'CholeskyPermutation', pcurr, ...
        'Display', false);
    
    % Typically the user wants the cholesky factorization
    varargout = {Rcurr, Pcurr};
else
    % Solve the problem using a reduced order model if the proper
    [Pf, L, d, B] = CAeigs(Aff, Bff, ne, Rold, Pold, ...
        dB(dofsfree, dofsfree), Xold(dofsfree, :), s, orthotype, orthovecs);
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