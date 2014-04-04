function [] = main11()

% two camera observing the origin
% a third camera is placed
% contourslice of the measure is plotted in 2D
%  - cartesian: x-y coordinates
%  - polar: d-a coordinates

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath

cartesian_or_polar = false;


if cartesian_or_polar % x-y
    gx = -2000:100:2000;
    gy = 10:100:2000;
    gz = 0:1:1;

    [nx,ny,nz] = meshgrid(gx,gy,gz);

    % measure, it has to be maximized
    nW = arrayfun(@(nx,ny,nz) min(eig(calc_Ciw( ...
           {CreateCamera('oripos', [0;0;1], [-1000; 500; 0]), ...
            CreateCamera('oripos', [0;0;1], [-1000; -500; 0]), ...
            CreateCamera('oripos', [0;0;1], [nx;ny;nz])}, ...
              [0;0;0]))), nx, ny, nz);

    figure
    hold on
    contourslice(nx,ny,nz,nW,[],[],0)
    shading flat
    axis('equal');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on
    hold off
    saveas(1, 'figures/slice_eig_newcam_cartesian.eps', 'epsc')
else % d-a
    gd = 10:10:2000;
    ga = 0:pi/10:pi;
    gb = 0:1:1;

    [nd,na,nb] = meshgrid(gd,ga,gb);

    % measure, it has to be maximized
    nW = arrayfun(@(nd,na,nb) min(eig(calc_Ciw( ...
           {CreateCamera('oripos', [0;0;1], [-1000; 500; 0]), ...
            CreateCamera('oripos', [0;0;1], [-1000; -500; 0]), ...
            CreateCamera('oripos', [0;0;1], [nd*sin(na);nd*cos(na);0])}, ...
              [0;0;0]))), nd, na, nb);

    figure
    hold on
    contourslice(nd,na,nb,nW,[],[],0)
    shading flat
    %axis('equal');
    xlabel('d');
    ylabel('a');
    zlabel('b');
    grid on
    hold off
    saveas(1, 'figures/slice_eig_newcam_polar.eps', 'epsc')
end

