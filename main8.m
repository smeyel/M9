function [] = main8()

% calculate the measure in a volume
% using virtual cameras with parameters defined in ICSSE_2013 publication
% the measure is the smallest eigenvalue of the inverse cov. matr.
% it has to be maximized
% the figure is contourslice

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


gx = -1000:100:1000;
gy = 100:100:1000;
gz = -1000:200:1000;

[nx,ny,nz] = meshgrid(gx,gy,gz);


cams = {CreateCamera('oripos', eye(3), [-1000; 0; 0]), ...
        CreateCamera('oripos', eye(3), [1000; 0; 0])};


% measure, it has to be maximized
nW = arrayfun(@(nx,ny,nz) min(eig(calc_Ciw(cams, [nx;ny;nz]))), nx, ny, nz);

figure
hold on
contourslice(nx,ny,nz,nW,[],[],gz)
shading flat
DrawCamera(cams{1}, 'r');
DrawCamera(cams{2}, 'g');
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
hold off

