function [] = main3()

% file for trying the actual problem
% There is given a camera configuration.
% There is given constraints for the new camera placement.
% Finding the optimal solution with fmincon.
% Plot the max-line from the partial derivatives.


%% preparation

close all
clear
clc

myAddPath

%colormap([zeros(63,3) ; ones(1,3)]);
%warning('off', 'Octave:possible-matlab-short-circuit-operator');

global useFoV;
useFoV=false;


%% calc


displayArea = [0 150 -60 70];


cam(1) = CreateCamera(-pi/4, [10;10]);
cam(2) = CreateCamera(pi/4, [10;-10]);

gX = 50;
gY = 0;

gsCovRes = CalculateResultingCovariance(cam, [gX;gY]);


maxSlope = getMaxSlope(gsCovRes, cam(1).e_pixel, cam(1).f_pixel);
maxRayFromSlope = createRay([gX gY], [gX+1 gY+maxSlope]);


[nX,nY] = meshgrid(85:1:100, -50:1:50);

nsCovRes = arrayfun(@(nx,ny) CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(gX-nx, gY-ny),[nx;ny])], [gX;gY]), ...
    nX, nY);
nW = arrayfun(@(covres) det(covres.Ci), nsCovRes);


minX = min(min(nX));
maxX = max(max(nX));
minY = min(min(nY));
maxY = max(max(nY));
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) myfunc(cam, gX, gY, x), ...
    [minX;minY], ...
    [], [], ...
    [], [], ...
    [minX;minY], [maxX;maxY]);

nArea = [minX minY ; ...
         minX maxY ; ...
         maxX maxY ; ...
         maxX minY];

maxXFromSolver = x;
maxRayFromSolver = createRay([gX gY], maxXFromSolver');


%figure
fig_my = figure; clf;
hold on
contour(nX,nY,nW, 60);
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
drawPolygon(nArea, 'k')

my_2D_error_ellipse(100*gsCovRes.C, [gX;gY], 'conf', 0.95);

plot(maxXFromSolver(1), maxXFromSolver(2), 'r*')
drawRay(maxRayFromSolver);

drawRay(maxRayFromSlope);

hold off



function f = myfunc(cam, gX, gY, x)

nx = x(1);
ny = x(2);

nsCovRes = CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(gX-nx, gY-ny),[nx;ny])], [gX;gY]);
nW = det(nsCovRes.Ci);

f = -nW;


function slope = getMaxSlope(covres, e, f)
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
slope = sqrt((A-B)/(A+B));
