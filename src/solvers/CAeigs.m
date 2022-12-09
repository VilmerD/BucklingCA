function [X, L] = CAeigs(A, B, n, Rold, Pold, dB, Xold, s, otype, ...
    ovecs)
% CAeigs(A, B, n, Rold, Pold, Bold, Xold, s, otype, ovecs) finds the n 
% eigenparis of the generalized eigenvalue problem correspodning to the n 
% eigenvalues with the largest absolute value.
% 
% Rold is the cholesky factorization B - dB. Xold are the old eigenvectors
% and s is the number of basis vectors to use. If s is a vector of length 
% n it uses variable number of basis vectors depending on the eigenpair. 
% If s is a scalar it uses s basis vectors for each eigenpair.
% 
% CAeigs orthogonalizes the basis vectors depending on the value of otype.
% If otype is...
% 
% CURRENT:
%       The current approximation of the eigenvectors are used.
% 
% OLD: 
%       The old approixmation of the eigenvectors in  ovecs is used.
%
% NONE or otherwise:
%       CA does not do basis orthogonalization. (NOT RECCOMENDED)

% Preprocess s
if numel(s) == 1; s = s*ones(n, 1); else; s = s(1:n); end

X = zeros(size(A, 1), n);
L = zeros(n, n);
for k = 1:n
    % Choose if orthogonalization is to be used and which vectors to
    % orthogonalize with respect to
    switch lower(otype)
        case 'current'
            XO = X(:, 1:(k-1));
        case 'old'
            XO = ovecs(:, 1:(k - 1));
        otherwise
            XO = [];
    end
    
    % Generate basis vectors
    Vk = CAEEON(A, B, Rold, Pold, dB, Xold(:, k), s(k), XO);
    
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
end

% Reorthogonalize vectors for 'old' case
if strcmpi(otype, 'old')
    X = reorthonormalize(X, B);
    L = diag(dot(X, A*X));
end
end

function X = reorthonormalize(X, A)
% Reorthogonalizes the vectors in X with respect to A
x1 = X(:, 1);
X(:, 1) = x1/sqrt(x1'*A*x1);
for k = 2:size(X, 2)
    % Extract element
    vk = X(:, k-1);

    % Project
    X(:, k:end) = X(:, k:end) - ((vk'*A)*X(:, k:end)).*vk;
    
    % Normalize 
    X(:, k) = X(:, k)/sqrt(X(:, k)'*A*X(:, k));
end

end