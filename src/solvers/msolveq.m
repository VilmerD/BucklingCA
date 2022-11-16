function [x, b, varargout] = msolveq(A, b, bc, Rold, Pold, dA, s)
% [x, b, R] = MSOLVEQ(A, b, bc) solves the linear system of equations 
% A*x = b with boundary conditions bc. b are the reaction forces and Rold is
% the cholesky factorization of the free part of A.
%
% [x, b, R] = MSOLVEQ(A, b, bc, R, P) solves the linear system of equations 
% where R is the cholesky factorization of A with permutation matrix P.
%
% [x, b, U, T, R, V] = MSOLVEQ(A, b, bc, Rold, dA, s) approximates the 
% solution u to the linear system of equations using CA. dA are the changes 
% in A since the last cholesky facotorization R. s A-orthonormal basis 
% vectors are computed using CASBON

% If bc is empty just solve directly
if isempty(bc)
    [Rcurr, FLAG, Pcurr] = chol(A, 'matrix');
    if FLAG
        errorstruct.message = 'A is not positive definite. Maybe forgot bc?';
        error(errorstruct);
    end
    x = Pcurr*(Rcurr\(Rcurr'\(Pcurr'*b)));
    return
end

% Split up the system into free and prescribed nodes
ndof = size(A, 1);
nf = (1:ndof)';
np = bc(:, 1);
nf(np) = [];

xp = bc(:, 2);
bf = b(nf) - A(nf, np)*xp;

% Solving free system
if nargin <= 5
    % Solve the problem using cholesky factorization
    if nargin < 4
        % No cholesky factorization of Aff is given, compute it
        [Rcurr, ~, Pcurr] = chol(A(nf, nf), 'matrix');
    elseif nargin == 4
        % No permutation matrix is given, assume no permutation
        Rcurr = Rold;
        Pcurr = speye(size(Rcurr));
    else
        % Both cholesky and P is given, just rename
        Rcurr = Rold;
        Pcurr = Pold;
    end
    xf = Pcurr*mldivide(Rcurr, mldivide(Rcurr', Pcurr'*bf));
    varargout = {Rcurr, Pcurr};
else
    % Using CA
    [V, U, T, R] = CASBON(A(nf, nf), bf, Rold, Pold, dA(nf, nf), s);
    xf = V*(V'*bf);
    varargout = {U, T, R, V};
end

% Reassembling the solution
x = zeros(ndof, 1);
x(np) = xp;
x(nf) = xf;

% Prescribed forces and reaction forces
b = zeros(ndof, 1);
b(np) = A(np, nf)*xf + A(np, np)*xp;
b(nf) = bf;
end