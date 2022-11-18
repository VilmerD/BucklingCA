%% Compute energy of buckling mode
% Strain energy?
% Pe'*KLe*Pe
% must iterate over elements
EPS = U(edofMat)*B';                                % strain
SIG = EPS*D;                                        % stress for E=1
EG = zeros(nelem, nevals);
for el = 1:nelem
    tau = [SIG(el,1) SIG(el,3)  0         0;
           SIG(el,3) SIG(el,2)  0         0;
           0         0          SIG(el,1) SIG(el,3);
           0         0          SIG(el,3) SIG(el,2)];
    BNLTtau = BNL'*tau;
    BNLTtauBNL = BNLTtau*BNL;
    eR = Emax*xPhys(el,1)^pS*helem^2;
%     eR = 1;
    elmdof = edofMat(el, :);
    Pe = PHI(elmdof, :);
    EG(el, :) = dot(Pe, BNLTtauBNL*Pe*eR);
end
% Stress energy?
% Pe'*KGe*Pe
% can be vectorized
EL = zeros(nelem, nevals);
for k = 1:nevals
    PHIk = PHI(:, k);
    Pe = PHIk(edofMat);
    eR = Emin+(xPhys').^pE*(Emax-Emin);
%     eR = 1;
    EL(:, k) = dot(Pe', (KE*Pe').*eR);
end

%% Plot some
fig = figure();
fig.Position(3:4) = [sizex, sizey]/helem;
ax = axes(fig, 'Position', [0, 0, 1, 1]);
k = 1;
xx = EG(:, k);
xx = EL(:, k);
xx = EG(:, k)./EL(:, k);
xx = (xx - min(xx))/(max(xx) - min(xx));
xx(xfilled) = -1;
% xx = xPhys;
mimg = sparse(i_img, j_img, xx);
o = imagesc(ax, mimg);
colormap(ax, 'hot');
%% Plot using fill
fig = figure();
fig.Position(3:4) = [sizex, sizey]/helem;
ax = axes(fig, 'Position', [0, 0, 1, 1]);
XX = X(:, 1)/helem;   YY = X(:, 2)/helem;
edofs = edofMat(:, 2:2:end)/2;
ex = XX(edofs); ey = YY(edofs);
xfilled = xPhys > 0.95;
exf = ex(xfilled, :);
eyf = ey(xfilled, :);
xxf = xx(xfilled);
fill(ax, exf', eyf', xxf, 'LineStyle', 'none');