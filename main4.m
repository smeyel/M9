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

calc_sym_s1_s2_K = true;
calc_x_y_not_a_d = false;

global firstCall


%% calculate the objective function

if calc_sym_s1_s2_K
syms s1 s2 s3 K4 K5 real
else
% actual values get from main3: gsCovRes.Ci
s1 = 20.4800;
s2 = 327.6800;
e = 1;
f = 480;
K = (e/f)^(-2);
end

if calc_x_y_not_a_d
    syms x y z real
    d = sqrt(x^2+y^2+z^2);
    a = atan(y/x);
    b = asin(z/d);
else
    syms a b d real
end

s4 = K4 * 1/d^2;
s5 = K5 * 1/d^2;
S1 = [s1 0 0 ; 0 s2 0 ; 0 0 s3];
S2 = [0 0 0 ; 0 s4 0 ; 0 0 s5];
Ra = [cos(a) -sin(a) 0 ; 0 0 1 ; sin(a) cos(a) 0];
Rb = [cos(b) 0 sin(b) ; 0 1 0 ; -sin(b) 0 cos(b)];
Se = S1 + Rb'*Ra'*S2*Ra*Rb

td = eig(Se)
d1 = td(1);
d2 = td(2);
d3 = td(3);

return

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


