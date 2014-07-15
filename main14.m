function [] = main14()

% fmincon on quasi-concave domains

close all
clear
clc


syms g1 g2 d1 d2 real
f = sin(g1-g2)^2 / (d1^2 * d2^2);

eps = 0.001;

count = 2;

%x = [g1;g2;d1;d2]
A = [0,0,-1,0; ...
    0,0,0,-1;...
    1,0,0,0;...
    -1,0,0,0;...
    0,1,0,0;...
    0,-1,0,0;...
    1,-1,0,0;...
    -1,1,0,0];
b = [-eps;-eps;pi/2;0;0;pi/2;pi;0];

[xopt q] = fmincon(...
        @(x) -myfun(reshape(x,length(x)/2,2)), ... %fun
        [zeros(count,1) ; 10*ones(count,1)], ... %x0
        A, b, ... %A, b
        [], [], ... %Aeq, beq
        [], ... %lb
        [], ... %ub
        @mycon, ... %nonlcon
        optimset('OutputFcn', @outputfun, 'MaxFunEvals', 3000)); %options
    %, 'Algorithm', 'interior-point'
xopt
q

return


function [c,ceq] = mycon(x)

m = 1;
a1 = 24*pi/180;
a2 = 20*pi/180;


g1 = x(1);
g2 = -x(2);
d1 = x(3);
d2 = x(4);
c(1) = m*sin(a1)/sin(a1+g1) - d1;
c(2) = m*sin(a2)/sin(a2+g2) - d2;
c
ceq = [];


function W = myfun(pols)
pol_cell = num2cell(pols,2);
[x y] = cellfun(@(pol) pol2cart(pol(1),pol(2)), pol_cell, 'uni', false);
cams = cellfun(@(x,y) CreateCamera('pos', [x;y]), x, y , 'uni', false);
Ci = calc_Ciw(cams, zeros(cams{1}.dim,1));
pols,Ci
%W = 1000*min(eig(Ci));
W = 1000*det(Ci);

g1 = pols(1,1);
g2 = pols(2,1);
d1 = pols(1,2);
d2 = pols(2,2);

g = [1/(d1^2*d2^2)*sin(2*(g1-g2));...
    1/(d1^2*d2^2)*sin(2*(g1-g2));...
    -2/(d1^3*d2^2)*sin(g1-g2)^2;...
    -2/(d1^2*d2^3)*sin(g1-g2)^2];
    

function stop = outputfun(t, optimValues, state)
stop = 0;


