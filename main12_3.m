function [] = main12_3()

% three cameras are placed on three parametric lines
% fmincon is applied
% steps and the resulting positions are plotted

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


global segment
segment{1} = [-500/cos(pi/3)+200, 0+200*tan(pi/6), 0, 500/cos(pi/6)];
segment{2} = [-500/cos(pi/3), 0, 500*cos(-2*pi/3),  500*sin(-2*pi/3)];
segment{3} = [500, 0, 700, 500];


figure
hold on
cellfun(@(s) plot([s(1);s(3)],[s(2);s(4)]), segment); % plot all segment
xopt = calc_opt_line_ori(segment);
plot(xopt(:,1), xopt(:,2), 'b*');
hold off


function [xopt q] = calc_opt_line_ori(segments)
x0 = cellfun(@(s) s(1), segments)';
y0 = cellfun(@(s) s(2), segments)';
x1 = cellfun(@(s) s(3), segments)';
y1 = cellfun(@(s) s(4), segments)';
dx = x1-x0;
dy = y1-y0;

[topt q] = fmincon(...
        @(t) -myfun([x0+t.*dx,y0+t.*dy]), ... %fun
        [0;0;0], ... %x0
        [], [], ... %A, b
        [], [], ... %Aeq, beq
        [0;0;0], ... %lb
        [1;1;1], ... %ub
        [], ... %nonlcon
        optimset('OutputFcn', @outputfun)); %options
xopt = [x0+topt.*dx,y0+topt.*dy];


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));


function stop = outputfun(t, optimValues, state)
% output function for the solver
% plots the camera positions

global segment
x0 = cellfun(@(s) s(1), segment)';
y0 = cellfun(@(s) s(2), segment)';
x1 = cellfun(@(s) s(3), segment)';
y1 = cellfun(@(s) s(4), segment)';
dx = x1-x0;
dy = y1-y0;

x = [x0+t.*dx,y0+t.*dy];

plot(x(1,1), x(1,2), 'g*')
plot(x(2,1), x(2,2), 'r*')
plot(x(3,1), x(3,2), 'm*')
stop = 0;


