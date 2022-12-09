%% Create a structured mesh for the bi-clamped beam
% X = nodal coordinates
% T = element connectivity and data
% row,col = for displaying design using imagesc (as a matrix)
function [X,T,row,col,solids,voids,F,freedofs] = generate_column(sizex,sizey,helem,doplot)
%% Define grid parameters 
dx = helem;                     % element size in x-direction
dy = helem;                     % element size in y-direction
loadL = 1/30*sizey;             % distribution length for point load
suppL = 1.00*loadL;
neumL = 1.00*suppL;
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
%% Boundary conditions for PDE filter
x_cent = elemlist_clean(:,1);
y_cent = elemlist_clean(:,2);
T(:,7:8) = [x_cent y_cent];
% Boundary have neumann
faceLB = x_cent-dx/2<100*eps;                   T(faceLB, 09) = 1;
faceBB = y_cent-dy/2<100*eps;                   T(faceBB, 10) = 1;
faceRB = x_cent+dx/2-sizex>-100*eps;            T(faceRB, 11) = 1;
faceTB = y_cent+dy/2-sizey>-100*eps;            T(faceTB, 12) = 1;
faceLB2 = bitand(x_cent-dx/2-sizex/2<100*eps, y_cent+dy/2-sizey/2>-100*eps);
faceTB2 = bitand(x_cent+dx/2-sizex/2<100*eps, y_cent+dy/2-sizey/2<100*eps);
T(faceLB2, 09) = 1;
T(faceTB2, 12) = 1;
% Support or load have robin
faceLS = bitand(...
    abs(x_cent-dx/2-0) < 100*eps, ...
    abs(y_cent-sizey/2)-neumL/2-dy/2<100*eps);  T(faceLS, 09) = 0;
faceRL = bitand(...
    abs(x_cent+dx/2-sizex) < 100*eps, ...
    abs(y_cent-sizey/2)-neumL/2-dy/2<100*eps);  T(faceRL, 11) = 0;
%% Define loads and supports 
solids = find(T(:,5)==2);   % find solid elements
voids = find(T(:,5)==1);    % find void elements

% Load at right end with x=sizex and y=sizey/2+-l_load
f0 = -1e4*helem;
loadednodes1 = bitand(...
    abs(X(:,1)-sizex)-0             <100*eps, ...
    abs(X(:,2)-sizey/2)-loadL/2     <100*eps);
loadednodes1 = find(loadednodes1);
loadeddofs1 = 2*loadednodes1-1;               % Force only in y
load1 = f0*ones(size(loadeddofs1));
loadednodes2 = bitand(....
    abs(X(:,1)-sizex)-0               <100*eps, ...
    bitor(...
        abs(X(:,2)-sizey/2-loadL/2-dy)-dy/2  <100*eps, ...
        abs(X(:,2)-sizey/2+loadL/2+dy)-dy/2  <100*eps));
loadednodes2 = find(loadednodes2);
loadeddofs2 = 2*loadednodes2-1;
load2 = (f0/2)*ones(size(loadeddofs2));
loadeddofs = [loadeddofs1; loadeddofs2];
load = [load1; load2];
F = sparse(loadeddofs,1,load,2*size(nodelist_clean,1),1,numel(loadeddofs));

% Supports at left end with x=0 and y=sizey/2+-l_load
supnodes1 = bitand(...
    abs(X(:,1)-0)-0                 <100*eps, ...
    abs(X(:,2)-sizey/2)-suppL/2     <100*eps);
supnodes1 = find(supnodes1);
supdofs1 = [2*supnodes1-1; 2*supnodes1];     % Support in y and y
% Supports at right end with x=sizex and y=sizey/2+-l_load
supnodes2 = bitand(...
    abs(X(:,1)-sizex)-0             <100*eps, ...
    abs(X(:,2)-sizey/2)-loadL/2     <100*eps);
supnodes2 = find(supnodes2);
supdofs2 = 2*supnodes2;     % Support only y
supdofs = unique([supdofs1;supdofs2]);

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

% Free dofs
alldofs = 1:2*size(nodelist_clean,1);
freedofs = setdiff(alldofs,supdofs)';
if (doplot); plot_geometry; end; end