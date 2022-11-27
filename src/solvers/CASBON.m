function V = CASBON(A, b, Rold, Pold, dA, s)
% V = CASBON(A, b, Rold, Pold, Aold, s) computes s orthonormal basis vectors using CA for
% the static problem A*x = b. R is the cholesky factorization of the old matrix Aold.
% 
% [V, U, T, R] = CASBON(A, b, Rold, Aold, s) computes the basis vectors and returns the 
% intermediate basis vectors required for a consistent sensitivity analysis.
n = size(b, 1);
AV = zeros(n, s);

% Compute first basis vector 
ui = Pold*mldivide(Rold, mldivide(Rold', Pold'*b));
ti = ui/sqrt(ui'*A*ui);
V = zeros(n, s);
V(:, 1) = ti;
AV(:, 1) = A*ti;

% Compute remaining basis vectors
for i = 2:s
    ui = -Pold*mldivide(Rold, mldivide(Rold', Pold'*(dA*ti)));
    ti = ui/sqrt(ui'*A*ui);
    
    % Orthogonalize
    ri = ti;
    for j = 1:(i-1)
        ri = ri - (ri'*AV(:, j))*V(:, j);
    end
    
    % Normalize
    vi = ri/sqrt(ri'*A*ri);
    
    % Insert
    V(:, i) = vi;
    AV(:, i) = A*vi;
end
end