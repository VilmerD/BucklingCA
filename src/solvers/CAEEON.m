function V = CAEEON(A, B, Rold, Pold, dB, Xold, s, X)
% V = CAEEON(A, B, Rold, Bold, Xold, s, XO) computes s basis vectors using CA for
% the eigenvalue problem (A - l*B)x = 0. R is the cholesky factorization of the 
% old matrix Bold. Xold are the old eigenvectors and X0 are vectors to which
% the basis vectors are orthogonalized. 
%
[n, m] = size(X);
BX = B*X;
BX = BX./dot(X, BX, 1);

% Compute first basis vector.
ui = Pold*mldivide(Rold, mldivide(Rold', Pold'*(A*Xold)));
ti = ui/sqrt(ui'*B*ui);
vi = ti;
for j = 1:m
    vi = vi - (ti'*BX(:, j))*X(:, j);
end
V = zeros(n, s);
V(:, 1) = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -Pold*mldivide(Rold, mldivide(Rold', Pold'*(dB*ti)));
    
    % Normalize
    ti = ui/sqrt(ui'*B*ui);
    
    % Orthogonalize
    vi = ti;
    for j = 1:m
        vi = vi - (ti'*BX(:, j))*X(:, j);
    end
    
    V(:, i) = vi;
end
end

%%% NOTES %%%
% - It's very important that X'*BX is I so that the orthogonalization is
% correct.
% - We use Grahm-Schmidt (not modified) since X may not be B-orthogonal