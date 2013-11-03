function [] = main4()

% file for trying the actual problem
% objective function: the smallest eigenvalue
%   of the inverse of the resulting covariance matrix
%   has to be maximized
% call the solver and plot the points


%% preparation

close all
clear
clc

myAddPath


%% definition

displayArea = [0 pi 0 100];

calc_sym = false;

if calc_sym
syms s1 s2 K
else
% actual values get from main3: gsCovRes.Ci
s1 = 20.4800;
s2 = 327.6800;
e = 1;
f = 480;
K = (e/f)^(-2);
end


%% calculate the objective function

syms a d

s3 = K * 1/d^2;
S1 = [s1 0 ; 0 s2];
S2 = [0 0 ; 0 s3];
R = [cos(a) -sin(a) ; sin(a) cos(a)];
Rt = [cos(a) sin(a) ; -sin(a) cos(a)];
Se = S1 + R*S2*Rt;

[tv,td] = eig(Se);
d1 = td(1,1);
d2 = td(2,2);

if calc_sym
obj = subs(d2, {K,s1,s2}, {(1/480)^(-2), 20.4800, 327.6800});
else
obj = d1;
end


%% base figures

fig_surf_eig = figure; clf;
ezsurf(obj, displayArea)
title('')

fig_contour_eig = figure; clf;
ezcontour(obj, displayArea)
title('')

fig_eig_fix_a_90deg = figure; clf;
ezplot(subs(obj, a, pi/2), [0 100])
title('')


%% Call the solver

minA = 0.4;
maxA = 2.6;
minD = 20;
maxD = 80;
startA = 0.5;
startD = 80;


oArea = [minA minD ; ...
         minA maxD ; ...
         maxA maxD ; ...
         maxA minD];


figure(fig_contour_eig);
hold on

drawPolygon(oArea)

global firstCall
firstCall = true;
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(...
    @(x) -subs(obj, {a,d}, {x(1),x(2)}), ... %fun
    [startA;startD], ... %x0
    [], [], ... %A, b
    [], [], ... %Aeq, beq
    [minA;minD], ... %lb
    [maxA;maxD], ... %ub
    [], ... %nonlcon
    optimset('OutputFcn', @outputfun)); %options

plot(x(1), x(2), 'r*')

hold off


%% save figures
saveas(fig_surf_eig, 'figures/surf_eig.eps', 'epsc')
saveas(fig_contour_eig, 'figures/contour_eig.eps', 'epsc')
saveas(fig_eig_fix_a_90deg, 'figures/eig_fix_a_90deg.eps', 'epsc')



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


