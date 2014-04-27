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


[px py] = pol2cart((-2:2:2)'*pi/3, 500);
[nx ny] = pol2cart((-2:2:2)'*pi/3+pi/2, 1);
d = [ -100 200 ; ...
      -300 600 ;...
      -400 400 ];
global segment
segment{1} = [px(1)+d(1,1)*nx(1), py(1)+d(1,1)*ny(1), px(1)+d(1,2)*nx(1), py(1)+d(1,2)*ny(1)];
segment{2} = [px(2)+d(2,1)*nx(2), py(2)+d(2,1)*ny(2), px(2)+d(2,2)*nx(2), py(2)+d(2,2)*ny(2)];
segment{3} = [px(3)+d(3,1)*nx(3), py(3)+d(3,1)*ny(3), px(3)+d(3,2)*nx(3), py(3)+d(3,2)*ny(3)];


global xopts
global qs
xopts = {};
qs = [];

calc_opt_line_ori(segment);

figure(1)
hold on
cellfun(@(s) plot([s(1);s(3)],[s(2);s(4)]), segment); % plot all segment
cellfun(@(x) plot(x(1,1), x(1,2), 'g*'), xopts) % plot cam 1 iteration positions
cellfun(@(x) plot(x(2,1), x(2,2), 'r*'), xopts) % plot cam 2 iteration positions
cellfun(@(x) plot(x(3,1), x(3,2), 'm*'), xopts) % plot cam 3 iteration positions
plot(xopts{end}(:,1), xopts{end}(:,2), 'b*'); %plot the optimal positions of all cameras
axis('equal')
hold off

figure(2)
plot(qs, 'k*')


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
        optimset('OutputFcn', @outputfun, 'MaxFunEvals', 300)); %options
xopt = [x0+topt.*dx,y0+topt.*dy];


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));


function stop = outputfun(t, optimValues, state)
% output function for the solver
% plots the camera positions

global segment
global xopts
global qs

if strcmp(state, 'iter')
    x0 = cellfun(@(s) s(1), segment)';
    y0 = cellfun(@(s) s(2), segment)';
    x1 = cellfun(@(s) s(3), segment)';
    y1 = cellfun(@(s) s(4), segment)';
    dx = x1-x0;
    dy = y1-y0;

    index = optimValues.iteration+1;
    x = [x0+t.*dx,y0+t.*dy];
    xopts{index} = x;
    qs(index) = myfun(x);
end
stop = 0;


