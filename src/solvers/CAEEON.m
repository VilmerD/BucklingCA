function [V, U, T] = CAEEON(A, B, Rold, Pold, Bold, Xold, s, XO)
% V = CAEEON(A, B, Rold, Bold, Xold, s, XO) computes s basis vectors using CA for
% the eigenvalue problem (A - l*B)x = 0. R is the cholesky factorization of the 
% old matrix Bold. Xold are the old eigenvectors and X0 are vectors to which
% the basis vectors are orthogonalized. 
% 
% [V, U, T] = CAEEON(A, B, Rold, Bold, Xold, s, XO) computes the basis vectors and returns the 
% intermediate basis vectors required for a consistent sensitivity analysis.
m = size(XO, 2);
dB = B - Bold;

% Compute first basis vector.
ui = Pold*(Rold\(Rold'\(Pold'*(A*Xold))));
ti = ui/sqrt(ui'*B*ui);
vi = ti;
for j = 1:m
    vi = vi - (vi'*A*XO(:, j))*XO(:, j);
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -Pold*(Rold\(Rold'\(Pold'*(dB*ti))));
    
    % Normalize
    ti = ui/sqrt(ui'*B*ui);
    
    % Orthogonalize
    vi = ti;
    for j = 1:m
        vi = vi - (vi'*A*XO(:, j))*XO(:, j);
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end
end