function [] = main2()

% file for trying the actual problem
% trying the solver call


%% notes
% single letters in variable names:
% n - new


%% preparation

close all
clear
clc

myAddPath

%colormap([zeros(63,3) ; ones(1,3)]);
%warning('off', 'Octave:possible-matlab-short-circuit-operator');

global useFoV;
useFoV=false;


%% local variable definitions

[nX,nY] = meshgrid(-49:2:49, -49:2:49);


%% Wellness with the new camera
 
nW = (-1) * arrayfun(@(nx,ny) myfunc([nx;ny]), nX, nY);

%figure
fig_contour_add_one_camera = figure; clf;
contour(nX,nY,nW,900:10:1100);
axis('equal');
xlabel('x');
ylabel('y', 'rotation', 0)


%% Call the solver

minX = -40;
maxX = -15;
minY = -30;
maxY = 30;
startX = -20;
startY = -20;

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@myfunc, ...
    [startX;startY], ...
    [], [], ...
    [], [], ...
    [minX;minY], [maxX;maxY]);

oArea = [minX minY ; ...
         minX maxY ; ...
         maxX maxY ; ...
         maxX minY];

%figure
figure(fig_contour_add_one_camera)
hold on
drawPolygon(oArea)
plot(x(1), x(2), 'r*')
hold off


%% save figures
saveas(fig_contour_add_one_camera, 'figures/contour_add_one_camera.eps')


function f = myfunc(x)
% There is given an origo centered covariance ellipse.
% There is given the reduced variance of the new camera.
% myfunc calculates the (-1)*Wellness of the new camera placement at
% the given position (x)

E = 10;
F = 90;
Gr = 0;
Hr = 1000;

t2 = x(1)^2 + x(2)^2;
K4 = (E-F)*(Gr-Hr) + Gr*Hr;
K2 = E*Hr + F*Gr;
K0 = E*F;

nW = x(2)^2 / t2^2 * K4 + ...
     1 / t2 * K2 + ...
     K0;
f = -nW;
