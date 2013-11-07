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

calc_sym_s1_s2_K = false;
calc_x_y_not_a_d = true;

global firstCall


%% calculate the objective function

if calc_sym_s1_s2_K
syms s1 s2 K real
else
% actual values get from main3: gsCovRes.Ci
s1 = 20.4800;
s2 = 327.6800;
e = 1;
f = 480;
K = (e/f)^(-2);
end

if calc_x_y_not_a_d
    syms x y real
    a = atan(y/x);
    d = sqrt(x^2+y^2);
else
    syms a d real
end

s3 = K * 1/d^2;
S1 = [s1 0 ; 0 s2];
S2 = [0 0 ; 0 s3];
R = [cos(a) -sin(a) ; sin(a) cos(a)];
Se = S1 + R*S2*R';

[tv,td] = eig(Se);
d1 = td(1,1);
d2 = td(2,2);

if calc_sym_s1_s2_K
obj = subs(d2, {K,s1,s2}, {(1/480)^(-2), 20.4800, 327.6800});
else
obj = d1;
end


%% base figures


if calc_x_y_not_a_d
    displayArea = [-100 100 0 100];

    fig_surf_eig_x_y = figure; clf;
    ezsurf(obj, displayArea)
    title('')

    fig_contour_eig_x_y = figure; clf;
    ezcontour(obj, displayArea)
    title('')

    fig_eig_fix_x_0 = figure; clf;
    ezplot(subs(obj, x, 0), [0 100])
    title('')
else
    displayArea = [0 pi 0 100];

    fig_surf_eig_a_d = figure; clf;
    ezsurf(obj, displayArea)
    title('')

    fig_contour_eig_a_d = figure; clf;
    ezcontour(obj, displayArea)
    title('')

    fig_eig_fix_a_90deg = figure; clf;
    ezplot(subs(obj, a, pi/2), [0 100])
    title('')
end


%% Call the solver (x, y)

if calc_x_y_not_a_d

minX = -50;
maxX = 50;
minY = 0;
maxY = 80;
startX = -40;
startY = 80;


oArea = [minX minY ; ...
         minX maxY ; ...
         maxX maxY ; ...
         maxX minY];


figure(fig_contour_eig_x_y);
hold on

drawPolygon(oArea)

firstCall = true;
xopt = fmincon(...
    @(v) -subs(obj, {x,y}, {v(1),v(2)}), ... %fun
    [startX;startY], ... %x0
    [], [], ... %A, b
    [], [], ... %Aeq, beq
    [minX;minY], ... %lb
    [maxX;maxY], ... %ub
    [], ... %nonlcon
    optimset('OutputFcn', @outputfun)); %options

plot(xopt(1), xopt(2), 'r*')

hold off

end


%% Call the solver (a, d)

if ~calc_x_y_not_a_d

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


figure(fig_contour_eig_a_d);
hold on

drawPolygon(oArea)

firstCall = true;
xopt = fmincon(...
    @(v) -subs(obj, {a,d}, {v(1),v(2)}), ... %fun
    [startA;startD], ... %x0
    [], [], ... %A, b
    [], [], ... %Aeq, beq
    [minA;minD], ... %lb
    [maxA;maxD], ... %ub
    [], ... %nonlcon
    optimset('OutputFcn', @outputfun)); %options

plot(xopt(1), xopt(2), 'r*')

hold off

end


%% save figures
if calc_x_y_not_a_d
    saveas(fig_surf_eig_x_y, 'figures/surf_eig_x_y.eps', 'epsc')
    saveas(fig_contour_eig_x_y, 'figures/contour_eig_x_y.eps', 'epsc')
    saveas(fig_eig_fix_x_0, 'figures/eig_fix_x_0.eps', 'epsc')
else
    saveas(fig_surf_eig_a_d, 'figures/surf_eig_a_d.eps', 'epsc')
    saveas(fig_contour_eig_a_d, 'figures/contour_eig_a_d.eps', 'epsc')
    saveas(fig_eig_fix_a_90deg, 'figures/eig_fix_a_90deg.eps', 'epsc')
end



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


