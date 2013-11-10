function [] = main5()

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
xlabel('x')
ylabel('z', 'rotation', 0)


%% save figures
saveas(fig_contour_3d, 'figures/contour_3d_1.eps', 'epsc')


function W = myfunc(v1, v2)

% parametric plane
p = [0;4;0];
p1 = [1;0;0];
p2 = [0;0;1];
xyz = p + v1*p1 + v2*p2;
x = xyz(1);
y = xyz(2);
z = xyz(3);

% actual values (s1,s2,e,f) get from main3: gsCovRes.Ci
s1 = 20.4800;
s2 = 327.6800;
s3 = 1.4*s2;
e = 1;
f = 480;
K4 = (e/f)^(-2);
K5 = K4 / 2;

d = sqrt(x^2+y^2+z^2);
a = GetAlpha2D(x,y);
b = asin(z/d);

s4 = K4 * 1/d^2;
s5 = K5 * 1/d^2;
S1 = [s1 0 0 ; 0 s2 0 ; 0 0 s3];
S2 = [0 0 0 ; 0 s4 0 ; 0 0 s5];
Ra = [cos(a) -sin(a) 0 ; 0 0 1 ; sin(a) cos(a) 0];
Rb = [cos(b) 0 sin(b) ; 0 1 0 ; -sin(b) 0 cos(b)];
Se = S1 + Rb'*Ra'*S2*Ra*Rb;
W = min(eig(Se));



