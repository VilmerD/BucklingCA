%% Plot nodes as points
figure;
axes(gcf);
hold on;
axis equal
axis tight
plot(nodelist_clean(:,1),nodelist_clean(:,2),'o', 'Markersize', 1);
% Plot supports and loads in x and y
supnodesx = (supdofs(logical(mod(supdofs, 2)))+1)/2;
plot(X(supnodesx, 1), X(supnodesx, 2), 'b>', 'Markersize', 6);
supnodesy = supdofs(logical(mod(supdofs-1, 2)))/2;
plot(X(supnodesy, 1), X(supnodesy, 2), 'b^', 'Markersize', 6);
loadednodesx = (loadeddofs(logical(mod(loadeddofs, 2)))+1)/2;
plot(X(loadednodesx, 1), X(loadednodesx, 2), 'r>', 'Markersize', 4);
loadednodesy = loadeddofs(logical(mod(loadeddofs-1, 2)))/2;
plot(X(loadednodesy, 1), X(loadednodesy, 2), 'r^', 'Markersize', 4);

% Plot elements
figure;
axes(gcf);
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

% Plot boundary elements for filter BCs
figure;
axes(gcf);
hold on;
axis equal
axis tight
eNeu = logical(sum(T(:, 9:12)==0, 2));
plot(T(eNeu, 7), T(eNeu, 8), 'm+');
eRob = logical(sum(T(:, 9:12)==1, 2));
plot(T(eRob, 7), T(eRob, 8), 'ro');
eDir = logical(sum(T(:, 9:12)==2, 2));