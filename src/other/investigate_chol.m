%% Compare using permutation matrix and not
N = 20;

t_factorize = 0;
t_factorize_with_perm = 0;

for k = 1:N
    tic;
    Rnp = chol(K(freedofs, freedofs));
    t_factorize = t_factorize + toc/N;

    tic;
    [Rp, FLAG, P] = chol(K(freedofs, freedofs));
    t_factorize_with_perm = t_factorize_with_perm + toc/N;
end
fprintf('Time to factorize without permutation: %f\n', t_factorize);
fprintf('Time to factorize with permutation: %f\n', t_factorize_with_perm);

tnp = 0;
tp = 0;
for k = 1:N
    tic;
    xknp = Rnp\(Rnp'\F(freedofs));
    t = toc;
    tnp = tnp + t/N;

    tic;
    xkp = P*(Rp\(Rp'\(P'*F(freedofs))));
    t = toc;
    tp = tp + t/N;

    if norm(xkp - xknp)/norm(xkp) > numel(xkp)*eps; error("Incorrect solution"); end
end

fprintf('Time without permutation: %f\n', tnp);
fprintf('Time with permutation: %f\n', tp);

%% Checking difference between vector and matrix permutation variants
Kff = K(freedofs, freedofs);
ff = F(freedofs);
[R, ~, P] = chol(Kff, 'matrix');
x1 = P*(R\(R'\(P'*ff)));
[R, ~, p] = chol(Kff, 'vector');
[~, pt] = sort(p);
x2 = R\(R'\ff(p));
x2 = x2(pt);

norm(x2 - x1)
norm(P'*Kff*P - Kff(p, p), 'fro')
%% How long time does it take to sort p, ie reverse it
t_sort = 0;
N = 1000;
for k = 1:N
    tic;
    [~, pt] = sort(p);
    t_sort = t_sort + toc/N;
end
fprintf('Time to sort p: %f\n', t_sort)