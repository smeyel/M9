function [] = main3()

% file for trying the actual problem
% There is given a camera configuration.
% There is given constraints for the new camera placement. It is built from
% one-parametrized curves.
% The contour of the objective function of the new camera placement is
% plotted.
% The max locations are calculated on every curve segment and plotted by
% red stars.
% The global maximum on the curves are plotted with green star.
% The global optimum found by the solver is plotted with a black star.


%% preparation

close all
clear
clc

myAddPath

%colormap([zeros(63,3) ; ones(1,3)]);
%warning('off', 'Octave:possible-matlab-short-circuit-operator');

global useFoV;
useFoV=false;


%% more preparations


displayArea = [0 150 -60 70];


cam(1) = CreateCamera(-pi/4, [10;10]);
cam(2) = CreateCamera(pi/4, [10;-10]);

gX = 50;
gY = 0;

gsCovRes = CalculateResultingCovariance(cam, [gX;gY]);


%% curves as constraints

syms a

f{1} = [85;-50]         + a * [0 ; 10];
f{2} = subs(f{1}, a, 1) + a * [5 ; 10];
f{3} = subs(f{2}, a, 1) + a * [-5 ; 10];
f{4} = subs(f{3}, a, 1) + a * [0 ; 20];
f{5} = subs(f{4}, a, 1) + a * [10 ; 10];
f{6} = subs(f{5}, a, 1) + a * [0 ; 40];
f{7} = subs(f{6}, a, 1) + a * [5 ; 0];
f{8} = subs(f{7}, a, 1) + a * [0 ; -100];
f{9} = subs(f{8}, a, 1) + a * [-15 ; 0];


%% objective function values for the contour plot

[nX,nY] = meshgrid(85:1:100, -50:1:50);
nsCovRes = arrayfun(@(nx,ny) CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(gX-nx, gY-ny),[nx;ny])], [gX;gY]), ...
    nX, nY);
nW = arrayfun(@(covres) det(covres.Ci), nsCovRes);


%% call the solver on the grid

% min or max on the new grid
minXn = min(min(nX));
maxXn = max(max(nX));
minYn = min(min(nY));
maxYn = max(max(nY));
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) objfunc(cam, gX, gY, x), ...
    [minXn;minYn], ...
    [], [], ...
    [], [], ...
    [minXn;minYn], [maxXn;maxYn]);
maxXFromSolver = x(1);
maxYFromSolver = x(2);


%% maxLines from the observed point

maxSlope = getMaxSlope(gsCovRes, cam(1).e_pixel, cam(1).f_pixel);
maxRayFromSlopeP = createRay([gX gY], [gX+1 gY+maxSlope]);
maxRayFromSlopeN = createRay([gX gY], [gX+1 gY-maxSlope]);


%% the extremal points of the curve segments

[objfuncc xc yc] = arrayfun(@(func) getObjFunc(gsCovRes, cam(1).e_pixel, cam(1).f_pixel, func), f, 'UniformOutput', false);
[maxAs maxXs maxYs maxFs] = cellfun(@(f, x, y) getObjFuncMaxForLine(f, x, y), objfuncc, xc, yc);
objfuncs = reshape([objfuncc{:}], size(objfuncc));


%% the global optimum

[maxF maxF_index] = max(maxFs);
maxX = maxXs(maxF_index);
maxY = maxYs(maxF_index);


%% figure

fig_my = figure; clf;
hold on
contour(nX,nY,nW, 60);
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
my_2D_error_ellipse(100*gsCovRes.C, [gX;gY], 'conf', 0.95);

cellfun(@(func) plotCurve(func), f);
plot(maxXs, maxYs, 'r*');
plot(maxX, maxY, 'g*');

drawRay(maxRayFromSlopeP);
drawRay(maxRayFromSlopeN);


plot(maxXFromSolver, maxYFromSolver, 'k*')

hold off



function plotCurve(func)
% plot one symbolic curve
syms a
b = 0:0.1:1;
Xc = arrayfun(@(b) subs(func, a, b), b, 'UniformOutput', false);
Xs = cell2mat(Xc);
plot(Xs(1,:), Xs(2,:));


function plotObjFuncOnCurveAndMax(func, maxA, maxF)
% plot the objective function on a curve and its extremal point
syms a
figure
hold on
ezplot(func, [0 1])
title('')
plot(maxA, maxF, 'r*')
hold off


function [f x y] = getObjFunc(covres, e, f, func)
% calculate the objective function

[A B C] = getConsts(covres, e, f);

syms a

xy0 = subs(func, a, 0);
xy1 = subs(func, a, 1);
x0 = xy0(1);
y0 = xy0(2);
x1 = xy1(1);
y1 = xy1(2);
dx = x1-x0;
dy = y1-y0;

x = x0+a*dx;
y = y0+a*dy;

% r: reduced, in the coordinate system fitted to the observed point
xr = x - covres.mu(1);
yr = y - covres.mu(2);
d2 = xr^2 + yr^2;
f = A * yr^2/d2^2 + B * 1/d2 + C;


function [maxA maxX maxY maxF] = getObjFuncMaxForLine(f, x, y)
% calculate the extremal point of the objective function on one curve

syms a
F=diff(f);
sols = subs(solve(F==0));
sols01 = sols(sols==real(sols) & 0 < sols & sols < 1);

as = [0;sols01;1];
fs=subs(f, a, as);
[maxF maxA_index] = max(fs);
maxA = as(maxA_index);
maxX = subs(x, a, maxA);
maxY = subs(y, a, maxA);


function f = objfunc(cam, gX, gY, x)
% objective function, it has to be minimized
nx = x(1);
ny = x(2);
nsCovRes = CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(gX-nx, gY-ny),[nx;ny])], [gX;gY]);
nW = det(nsCovRes.Ci);
f = -nW;


function [A B C] = getConsts(covres, e, f)
% calculate the constant values A, B and C
% e, f: error and focal length of the new camera, same unit
[v,d] = eig(covres.Ci);
if(d(1,1)<d(2,2))
    s1 = d(1,1);
    s2 = d(2,2);
else
    s1 = d(2,2);
    s2 = d(1,1);
end
A = (s2-s1) * (e/f)^(-2);
B = s1 * (e/f)^(-2);
C = s1 * s2;


function slope = getMaxSlope(covres, e, f)
% calculate the slope of the maxLine from the observed point
[A B C] = getConsts(covres, e, f);
slope = sqrt((A-B)/(A+B));
