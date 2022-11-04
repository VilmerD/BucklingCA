%% Create a structured mesh for the bi-clamped beam
% X = nodal coordinates
% T = element connectivity and data
% row,col = for displaying design using imagesc (as a matrix)
function [X,T,row,col,solids,voids,F,freedofs] = generate_columnB(sizex,sizey,helem,doplot)
%% Define grid parameters 
dx = helem;                     % element size in x-direction
dy = helem;                     % element size in y-direction
l_load = 0.04*sizey;            % distribution length for point load
l_supp = 0;
l_neu = l_load+dy;
%% Create nodal grid for FEA and optimization
[Xnode,Ynode] = meshgrid(0:dx:sizex,0:dy:sizey);
nodelist(:,1) = Xnode(:);
nodelist(:,2) = Ynode(:);
xnodeout = [];
ynodeout = [];
xynodeout = intersect(xnodeout,ynodeout);
nodelist_clean = nodelist;
nodelist_clean(xynodeout,:) = [];
X = nodelist_clean; % For output
%% Create element grid for FEA and optimization
[Xelem,Yelem] = meshgrid(dx/2:dx:sizex-dx/2,dy/2:dy:sizey-dy/2);
elemlist(:,1) = Xelem(:);
elemlist(:,2) = Yelem(:);
xelemout = [];
yelemout = []; 
xyelemout = intersect(xelemout,yelemout);
elemlist_clean = elemlist;
elemlist_clean(xyelemout,:) = [];
% Create element connectivity
nelem = size(elemlist_clean,1);
T = zeros(nelem,12);
% T(1:4) = [node1 node2 node3 node4]
% T(5) = 0/1/2 (design / void / solid)
% T(6) = 0/1 (domain / padding region)
% T(7:8) = [x_centroid y_centroid]
% T(9:12) = [left bottom right top] BCs for PDE filter
% This is a rectangular mesh - easy to create based on X and Y
for ex = 1:size(Xelem,2)
    for ey = 1:size(Yelem,1)
        e = (ex-1)*size(Xelem,1)+ey;
        n1 = (ex-1)*size(Xnode,1)+ey;
        n2 = n1+size(Xnode,1);
        n3 = n2+1;
        n4 = n1+1;
        T(e,1:4) = [n1 n2 n3 n4];
    end
end
for e = 1:nelem
    % Centroid coordinates
    x_cent = elemlist_clean(e,1);
    y_cent = elemlist_clean(e,2);
    T(e,7:8) = [x_cent y_cent];
    % Assign voids
    % No voids in the current representation of the column
    % Assign solids (for loaded and supported regions)
%   if (x_cent>sizex-helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,5:6) = [2 0]; end
% 	if (x_cent<helem && y_cent>sizey/2-l_load/2 && y_cent<sizey/2+l_load/2); T(e,5:6) = [2 0]; end
    % Assign non-Neumann borders (for augmented PDE filter)
    % 0 = Neumann (supports, load)
    % 1 = Robin BCs with l_s=l_o (free faces)
    % 2 = xi_corner (Dirichlet-like)
    % Default is Neumann [left bottom right top]
    if (x_cent-dx/2<100*eps);          T(e,09) = 1; end    % left face
    if (y_cent-dy/2<100*eps);          T(e,10) = 1; end    % bottom face
    if (x_cent+dx/2-sizex>-100*eps);   T(e,11) = 1; end    % right face
    if (y_cent+dy/2-sizey>-100*eps);   T(e,12) = 1; end    % top face
    if (abs(x_cent-dx/2 - 0)     < 100*eps  && abs(y_cent-sizey/2)-l_neu/2<100*eps); T(e,09) = 0; end % left face, near support
    if (abs(x_cent+dx/2 - sizex) < 100*eps  && abs(y_cent-sizey/2)-l_neu/2<100*eps); T(e,11) = 0; end % right face, near load
end
%% Create matrix representation of topology
xmin = dx/2; 
ymin = dy/2; ymax = sizey-dy/2;
nrows = (ymax-ymin)/dy + 1;
col = zeros(nelem,1);
row = col;
for e = 1:nelem
    col(e) = (T(e,7)-xmin)/dx + 1;
    row(e) = nrows - ((T(e,8)-ymin)/dy);
end
%% Define loads and supports 
solids = find(T(:,5)==2);   % find solid elements
voids = find(T(:,5)==1);    % find void elements

% Load at right end with x=sizex and y=sizey/2+-l_load
loadednodes1 = bitand(...
    abs(X(:,1)-sizex)-0             <100*eps, ...
    abs(X(:,2)-sizey/2)-l_load/2    <100*eps);
loadednodes1 = find(loadednodes1);
loadeddofs1 = 2*loadednodes1-1;               % Force only in y
load1 = -1e6/numel(loadeddofs1)*ones(size(loadeddofs1));

% Load at right end with x=sizex and y=0+-l_load
loadednodes2 = bitand(...
    abs(X(:,1)-0)-0                 <100*eps, ...
    abs(X(:,2)-sizey/2)-l_load/2    <100*eps);  
loadednodes2 = find(loadednodes2);
loadeddofs2 = 2*loadednodes2-1;               % Force only in y
load2 = +1e6/numel(loadeddofs2)*ones(size(loadeddofs2));
loadeddofs = union(loadeddofs1, loadeddofs2);
load = [load1; load2];
F = sparse(loadeddofs,1,load,2*size(nodelist_clean,1),1,numel(loadeddofs));

% Supports at left end with x=0 and y=sizey/2+-l_load
supnodes1 = bitand(...
    abs(X(:,1)-0)-0                 <100*eps, ...
    abs(X(:,2)-sizey/2)-l_supp/2    <100*eps);
supnodes1 = find(supnodes1);
supdofs1 = 2*supnodes1;     % Support in y
% Supports at right end with x=sizex and y=sizey/2+-l_load
supnodes2 = bitand(...
    abs(X(:,1)-sizex)-0             <100*eps, ...
    abs(X(:,2)-sizey/2)-l_supp/2    <100*eps);
supnodes2 = find(supnodes2);
supdofs2 = 2*supnodes2;     % Support only y
% Support in the middle with x = sizex/2 and y = sizey/2
supnodes3 = bitand(...
    abs(X(:,1)-sizex/2)-0       <100*eps, ...
    abs(X(:,2)-sizey/2)-0       <100*eps);
supnodes3 = find(supnodes3);
supdofs3 = 2*supnodes3-1;   % Support only x
supdofs = unique([supdofs1;supdofs2;supdofs3]);

% Free dofs
alldofs = 1:2*size(nodelist_clean,1);
freedofs = setdiff(alldofs,supdofs)';
%% Plot nodes as points
if (doplot)
    figure;
    plot(nodelist_clean(:,1),nodelist_clean(:,2),'o');
    hold on;
    axis equal
    axis tight
    % Plot supports and loads in x and y
    supnodesx = (supdofs(logical(mod(supdofs, 2)))+1)/2;
    plot(X(supnodesx, 1), X(supnodesx, 2), 'b>');
    supnodesy = supdofs(logical(mod(supdofs-1, 2)))/2;
    plot(X(supnodesy, 1), X(supnodesy, 2), 'b^');
    loadednodesx = (loadeddofs(logical(mod(loadeddofs, 2)))+1)/2;
    plot(X(loadednodesx, 1), X(loadednodesx, 2), 'r>');
    loadednodesy = loadeddofs(logical(mod(loadeddofs-1, 2)))/2;
    plot(X(loadednodesy, 1), X(loadednodesy, 2), 'r^');
end
% Plot elements
if (doplot)
    figure;
    hold on;
    axis equal
    axis tight
    % Regular element
    eRegular = T(e, 5) == 0;
    plot(T(eRegular, 7), T(eRegular, 8), 'ro');
    % Void element
    eVoid = T(e, 5) == 1;
    plot(T(eVoid, 7), T(eVoid, 8), 'm+');
    % Solid element
    eSolid = T(e, 5) == 2;
    plot(T(eSolid, 7), T(eSolid, 8), 'k*');
end
% Plot boundary elements for filter BCs
if (doplot)
    figure;
    hold on;
    axis equal
    axis tight
    eNeu = logical(sum(T(:, 9:12)==0, 2));
    plot(T(eNeu, 7), T(eNeu, 8), 'm+');
    eRob = logical(sum(T(:, 9:12)==1, 2));
    plot(T(eRob, 7), T(eRob, 8), 'ro');
    eDir = logical(sum(T(:, 9:12)==2, 2));
    plot(T(eDir, 7), T(eDir, 8), 'k*');
end
end