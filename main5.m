function [] = main5()

% ----------
% OBSOLETE BECAUSE OF THE CAMERA COORD. SYSTEM CONVERSION
% FACING DIRECTION IS THE LAST DIMENSION, NOT THE X
% ----------
% file for trying the actual problem
% 3D
% objective function: the smallest eigenvalue
%   of the inverse of the resulting covariance matrix
%   has to be maximized
% calculation is done on a meshgrid by using arrayfun
% contour is plotted


%% preparation

close all
clear
clc

myAddPath


%% calculate the objective function
[n1,n2] = meshgrid(-39:39, -39:39);
nW = arrayfun(@(n11,n22) myfunc(n11, n22), n1, n2);


%% figure
fig_contour_3d = figure; clf;
contour(n1, n2, nW);
axis('equal');
%xlabel('z')
%ylabel('x', 'rotation', 0)


%% save figures
colormap([zeros(63,3) ; ones(1,3)]);
saveas(fig_contour_3d, 'figures/contour_3d_1.eps', 'epsc')


function W = myfunc(d1, d2)

% parametric plane
p = [0;4;0];
w1 = [0;0;1];
w2 = [1;2;0];
v1 = w1 / norm(w1);
w2n = w2 - v1'*w2*v1;
v2 = w2n / norm(w2n);
xyz = p + d1*v1 + d2*v2;
x = xyz(1);
y = xyz(2);
z = xyz(3);

% actual values (s1,s2,e,f) get from main3: gsCovRes.Ci
s1 = 20.4800;
s2 = 327.6800;
s3 = 2*s2;
e = 1;
f = 480;
K4 = (e/f)^(-2);
K5 = K4 / 2;

[a b d] = cart2sph(x, y, z);
b = -b; % elevation from xy plane, but clockwise around y axis
c = 1;

s4 = K4 * 1/d^2;
s5 = K5 * 1/d^2;
S1 = [s1 0 0 ; 0 s2 0 ; 0 0 s3];
S2 = [0 0 0 ; 0 s4 0 ; 0 0 s5];
Rx = [1 0 0 ; 0 cos(c) -sin(c) ; 0 sin(c) cos(c)];
Ry = [cos(b) 0 sin(b) ; 0 1 0 ; -sin(b) 0 cos(b)];
Rz = [cos(a) -sin(a) 0 ; sin(a) cos(a) 0 ; 0 0 1];
Se = S1 + Rz*Ry*Rx*S2*Rx'*Ry'*Rz';
W = min(eig(Se));



