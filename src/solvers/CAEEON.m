function [V, U, T] = CAEEON(A, B, Rold, Bold, Pold, s, VO)
% V = CAEEON(A, B, Rold, Bold, Pold, s, VO) computes s basis vectors using CA for
% the eigenvalue problem (A - l*B)x = 0. R is the cholesky factorization of the 
% old matrix Bold. Pold are the old eigenmodes and V0 are vectors to which
% the basis vectors are orthogonalized. 
% 
% [V, U, T] = CAEEON(A, B, Rold, Bold, Pold, s, VO) computes the basis vectors and returns the 
% intermediate basis vectors required for a consistent sensitivity analysis.
m = size(VO, 2);
dB = B - Bold;

% Compute first basis vector.
ui = Rold\(Rold'\(A*Pold));
ti = ui/sqrt(ui'*B*ui);
vi = ti;
for j = 1:m
    vj = VO(:, j);
    vi = vi - (ti'*A*vj)*vj;
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -Rold\(Rold'\(dB*ti));
    
    % Normalize
    ti = ui/sqrt(ui'*B*ui);
    
    % Orthogonalize
    vi = ti;
    for j = 1:m
        vj = VO(:, j);
        vi = vi - (ti'*A*vj)*vj;
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end
end