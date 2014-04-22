function [] = main12()

% two camera observing the origin
% a third camera is placed
% fmincon in 2D

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

global firstCall

myAddPath

cams = {CreateCamera('oripos', [0;1], [-1000; 500]), ...
        CreateCamera('oripos', [0;1], [-1000; -500])};
    
constraints = [ -500, 500 ; ...
                -500, 2000 ; ...
                 500, 2000 ; ...
                 500, 500 ];
constraints = [constraints ; constraints(1,:)];

F = (cams{1}.e/cams{1}.fx)^(-2);
ee =  eig(calc_Ciw(cams, zeros(cams{1}.dim,1)));
E = ee(2)-ee(1);
D = ee(2)+ee(1);
dist = sqrt(F/E);

            
hold on
DrawCamera(cams)
plot(constraints(:,1),constraints(:,2))

start = [500; 1500];
lb = [0;10];
ub = [5000;5000];

firstCall = true;
xopt = fmincon(...
    @(v) -myfun(cams,v), ... %fun
    start, ... %x0
    [], [], ... %A, b
    [], [], ... %Aeq, beq
    lb, ... %lb
    ub, ... %ub
    [], ... %nonlcon
    optimset('OutputFcn', @outputfun)); %options
plot(xopt(1), xopt(2), 'r*')
hold off

xopt
myfun(cams, xopt)


if inpolygon(xopt(1), xopt(2), constraints(:,1),constraints(:,2))
    disp('IN')
else
    disp('OUT')
end


return


function W = myfun(cams, x)
cams{length(cams)+1} = CreateCamera('oripos', [0;1], x);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));



function stop = outputfun(x, optimValues, state)
% output function for the solver
% plots the points and connects the consecutive ones
global firstCall
global pre
if ~firstCall
    plot([x(1) pre(1)], [x(2) pre(2)], 'k')
end
firstCall = false;
plot(x(1), x(2), 'g*')
pre = x;
stop = 0;

