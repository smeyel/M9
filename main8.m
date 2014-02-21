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
scatter3(cams{1}.pos(1), cams{1}.pos(2), cams{1}.pos(3), 'r', 'filled');
scatter3(cams{2}.pos(1), cams{2}.pos(2), cams{2}.pos(3), 'g', 'filled');
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
hold off

