function [] = main6()

% calculate the resulting covariance matrix
% in case of one specific observed point
% with the cameras and measurement setup
% in the ICSSE_2013 publication

close all
clear
clc

myAddPath

% camera positions
% from the ICSSE_2013 publication
% from the log
t0 = [-674.3543 ; 263.9349 ; -677.6962];
t1 = [871.0314 ; -266.2977 ; -522.8290];
t2 = [-28.0911 ; -265.4436 ; -758.5166];

% observed point
p = [175.09 ; -185.7 ; 103.74];

% covariance matrix inverses and the resulting covariace matrix
Ciw0 = calc_Ci(t0, p);
Ciw1 = calc_Ci(t1, p);
Ciw2 = calc_Ci(t2, p);
Ciw = Ciw0 + Ciw1 + Ciw2;
Cw = inv(Ciw);

% variance in the x-y-z directions
[V D] = eig(Cw);
sigres = V.^2 * diag(D)

% figure
figure
hold on
error_ellipse(Cw, p, 'conf', 0.95);
scatter3(t0(1),t0(2),t0(3), 'r', 'filled');
scatter3(t1(1),t1(2),t1(3), 'g', 'filled');
scatter3(t2(1),t2(2),t2(3), 'b', 'filled');
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
hold off


% calculate the inverse of the covariance matrix
% for one camera and one observed point
% the result is given in the world coordinate system
function Ciw = calc_Ci(t, p)

e = 1; % assumption
f = 789.1510; % from the ICSSE_2013 publication, from the C++ code

v = p-t;
d = norm(v);
sig = d * e / f;
si = sig^(-2);
Ci = diag([0,si,si]);

x = 0; % fx = fy = f, symmetric
y = asin(v(3)/d);
z = GetAlpha2D(v(1),v(2));

Rz = [cos(z) -sin(z) 0 ; sin(z) cos(z) 0 ; 0 0 1];
Ry = [cos(y) 0 sin(y) ; 0 1 0 ; -sin(y) 0 cos(y)];
Rx = [1 0 0 ; 0 cos(x) -sin(x) ; 0 sin(x) cos(x)];
Ciw = Rz*Ry*Rx*Ci*Rx'*Ry'*Rz';
