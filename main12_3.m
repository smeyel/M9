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

num = 4; % 3 or 4
[px py] = pol2cart((floor(-num/2+1):floor(num/2))'*(2*pi/num), 500);
[nx ny] = pol2cart((floor(-num/2+1):floor(num/2))'*(2*pi/num)+pi/2, 1);
d = [ -100 200 ; ...
      -300 600 ;...
      -400 400 ;...
      -600 600 ];
cc = {'g*', 'r*', 'm*', 'c*'};
global segment
segment = cell(1,num);
for n=1:num
    segment{n} = [px(n)+d(n,1)*nx(n), py(n)+d(n,1)*ny(n), px(n)+d(n,2)*nx(n), py(n)+d(n,2)*ny(n)];
end

segment{1} = [-500 -400 400 -600];
segment{2} = [600 -300 -200 400];
segment{3} = [-700 0 200 600];
segment{4} = [-400 -400 -400 400];


global xopts
global qs
global funccount
xopts = {};
qs = [];
funccount = [];

calc_opt_line_ori(segment);

figure(1)
hold on
cellfun(@(s) plot([s(1);s(3)],[s(2);s(4)]), segment); % plot all segment
for n=1:num
    cellfun(@(x) plot(x(n,1), x(n,2), cc{n}), xopts) % plot cam 'n' iteration positions with color cc{n}
end
plot(xopts{end}(:,1), xopts{end}(:,2), 'b*'); %plot the optimal positions of all cameras
axis('equal')
hold off

figure(2)
hold on
plot(qs)
plot(qs, 'k*')
hold off

figure(3)
hold on
plot(funccount)
plot(funccount, 'k*')
hold off


function [xopt q] = calc_opt_line_ori(segments)
x0 = cellfun(@(s) s(1), segments)';
y0 = cellfun(@(s) s(2), segments)';
x1 = cellfun(@(s) s(3), segments)';
y1 = cellfun(@(s) s(4), segments)';
dx = x1-x0;
dy = y1-y0;

count = numel(segments);

[topt q] = fmincon(...
        @(t) -myfun([x0+t.*dx,y0+t.*dy]), ... %fun
        zeros(count,1), ... %x0
        [], [], ... %A, b
        [], [], ... %Aeq, beq
        zeros(count,1), ... %lb
        ones(count,1), ... %ub
        [], ... %nonlcon
        optimset('OutputFcn', @outputfun, 'MaxFunEvals', 300)); %options
xopt = [x0+topt.*dx,y0+topt.*dy];


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));


function stop = outputfun(t, optimValues, state)
% output function for the solver
% plots the camera positions

global segment
global xopts
global qs
global funccount

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
    qs(index) = 1/myfun(x);
    funccount(index) = optimValues.funccount;
end
stop = 0;


