function [X, L, V_INTR] = CAeigs(A, B, n, Rold, Pold, dB, Xold, s, orthotype, ...
    orthovecs)
% CAeigs(A, B, n, Rold, Pold, Bold, Xold, s, options) finds the n eigenparis of the
% generalized eigenvalue problem correspodning to the n eigenvalues with
% the largest absolute value.
% 
% Rold is the cholesky factorization of Bold, an older version of B 
% Xold are the old eigenvectors and s is the number of basis vectors to use.
% 
% Finally options is a struct containing the orthotype and orthovecs
% 
% If orthotype is CURRENT CA uses the current approximation of the
% eigenvectors for orthogonalization.
% 
% If orthotype is OLD CA uses an old approixmation of the eigenvectors
% for orthogonalization, which must be contained in orthovecs
%
% If orthotype is NONE or otherwise, CA does not do basis
% orthogonalization.

% Preprocess s
if numel(s) == 1; s = s*ones(n, 1); else; s = s(1:n); end

X = zeros(size(A, 1), n);
L = zeros(n, n);
V_INTR = cell(n, 4);
for k = 1:n
    % Choose if orthogonalization is to be used and which vectors to
    % orthogonalize with respect to
    switch lower(orthotype)
        case 'current'
            XO = X(:, 1:(k-1));
        case 'old'
            XO = orthovecs(:, 1:(k - 1));
        otherwise
            XO = [];
    end
    
    % Generate basis vectors
    [Vk, Uk, Tk] = CAEEON(A, B, Rold, Pold, dB, Xold(:, k), s(k), XO);
    
    % Compute reduced model
    A_RED = Vk'*A*Vk;
    B_RED = Vk'*B*Vk;
    
    % Solve reduced problem
    [X_RED, L_RED] = eigs(A_RED, B_RED, 1, ...
        'largestreal', ...
        'IsCholesky', false);
    
    % Normalize the eigenvectors with respect to B
    X_RED = X_RED/sqrt(X_RED'*B_RED*X_RED);                    
    X_FULL = Vk*X_RED;
    
    % Insert solution
    X(:, k) = X_FULL;                             
    L(k, k) = L_RED;
    V_INTR(k, :) = {Uk, Tk, Vk, X_RED};
    
    % Compute the residual of the eigenproblem
%     d(k) = norm(A*X_FULL - L_RED*B*X_FULL)/norm(A*X_FULL);
end

end