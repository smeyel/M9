function [] = main12_2()

% one camera observing the object at origin
% a second and third camera are placed on two parametric lines
% contour is plotted

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


s1 = [-500/cos(pi/3)+200, 0+200*tan(pi/6), 0, 500/cos(pi/6)];
s2 = [-500/cos(pi/3), 0, 500*cos(-2*pi/3),  500*sin(-2*pi/3)];
c1 = CreateCamera('pos', [500; 0]);


g1 = linspace(0,1,50);
g2 = linspace(0,1,50);
[t1,t2] = meshgrid(g1,g2);

x01 = s1(1);
y01 = s1(2);
x11 = s1(3);
y11 = s1(4);
dx1 = x11-x01;
dy1 = y11-y01;

x02 = s2(1);
y02 = s2(2);
x12 = s2(3);
y12 = s2(4);
dx2 = x12-x02;
dy2 = y12-y02;

tW = arrayfun(@(t1,t2) min(eig(calc_Ciw( ...
           {c1, ...
            CreateCamera('pos', [x01+t1*dx1;y01+t1*dy1]), ...
            CreateCamera('pos', [x02+t2*dx2;y02+t2*dy2])}, ...
              [0;0]))), t1, t2);

figure
hold on
surf(t1,t2,tW)
xlabel('t1');
ylabel('t2');
grid on
hold off


