function [] = main12_5()

% 2D parallel optimization with fmincon with half plane constraints

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


num = 2;
cc = {'g*', 'r*', 'm*', 'c*'};


global xopts
global qs
xopts = {};
qs = [];

A = [-1  0  0  0 ; ...
      0  0  0 -1 ; ...
      0  0 -1  1 ];
b = [ -500 ; ...
      -500 ; ...
      1000 ];

[xopt q] = fmincon(...
        @(xa) -myfun([xa(1:2:end),xa(2:2:end)]), ... %fun
        [600;-700;700;900], ... %x0
        A, b, ... %A, b
        [], [], ... %Aeq, beq
        [], ... %lb
        [], ... %ub
        [], ... %nonlcon
        optimset('OutputFcn', @outputfun, 'MaxFunEvals', 300)); %options




figure(1)
hold on
for n=1:num
    cellfun(@(x) plot(x(n,1), x(n,2), cc{n}), xopts) % plot cam 'n' iteration positions with color cc{n}
end
plot(xopts{end}(:,1), xopts{end}(:,2), 'b*'); %plot the optimal positions of all cameras
axis('equal')
hold off

figure(2)
plot(qs, 'k*')

return




function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));


function stop = outputfun(xa, optimValues, state)
% output function for the solver
% plots the camera positions

global xopts
global qs

if strcmp(state, 'iter')
    index = optimValues.iteration+1;
    x = [xa(1:2:end),xa(2:2:end)];
    xopts{index} = x;
    qs(index) = myfun(x);
end
stop = 0;


