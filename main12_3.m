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

num = 2;
segment = cell(1,num);
%segment{3} = [-700 0 200 600];
%segment{4} = [-400 -400 -400 400];
%segment{1} = [-500 -270 400 -270];
%segment{2} = [600 -300 -200 500];
k = -150*sqrt(2);
k = -300;
segment{1} = [0 k -k k];
segment{2} = [150 150 300 0];

o1 = [0 k ; 300 0];
o2 = [k -k ; 150 150];


m = 1;
num = 2;
segment = cell(1,num);
a1 = 60*pi/180;
a2 = 80*pi/180;
segment{1} = [m 0 0 m*tan(a1)];
segment{2} = [m 0 0 -m*tan(a2)];



global xopts
global qs
global funccount
xopts = {};
qs = [];
funccount = [];

[xopt q] = calc_opt_line_ori(segment)

figure(1)
hold on
cellfun(@(s) plot([s(1);s(3)],[s(2);s(4)]), segment); % plot all segment
for n=1:num
    cellfun(@(x) plot(x(n,1), x(n,2), cc{n}), xopts) % plot cam 'n' iteration positions with color cc{n}
end
plot(xopts{end}(:,1), xopts{end}(:,2), 'b*'); %plot the optimal positions of all cameras
plot(0, 0, 'k*')
axis('equal')
hold off

%myfun(o1)
%myfun(o2)
x1 = xopt(1,:)';
x2 = xopt(2,:)';
x1'*x2
norm(x1)
norm(x2)
return

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
topt
xopt = [x0+topt.*dx,y0+topt.*dy];


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
Ci = calc_Ciw(cams, zeros(cams{1}.dim,1));
W = min(eig(Ci));
% W = -trace(pinv(Ci));
% W = det(Ci);


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


