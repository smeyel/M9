function [] = main7()

% calculate the measure
% in the volume defined in the ICSSE_2013 publication
% the measure is the smallest eigenvalue of the inverse cov. matr.
% it has to be maximized
% the cameras and measurement setup is defined in ICSSE_2013 publication
% the figure is slice

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath

indata = dlmread('M1R1_stableframeDataProcessor_stats.csv',' ');

% Statistics for every stable frame interval (location)
LocationMean2Ray = indata(:,1:3);
LocationMean3Ray = indata(:,4:6);
LocationMeanAll = indata(:,7:9);
LocationStd2Ray = indata(:,10:12);
LocationStd3Ray = indata(:,13:15);
LocationStdAll = indata(:,16:18);
LocationEffectiveStd2Ray = indata(:,19);
LocationEffectiveStd3Ray = indata(:,20);
LocationEffectiveStdAll = indata(:,21);


minLoc = min(LocationMeanAll)';
maxLoc = max(LocationMeanAll)';

step = 10;
gx = minLoc(1):step:maxLoc(1);
gy = minLoc(2):step:maxLoc(2);
gz = minLoc(3):step:maxLoc(3);
[nx,ny,nz] = meshgrid(gx,gy,gz);


% from the ICSSE_2013 publication
% T matrices from the cpp program print
% transformation: camera => world
m0 = [-0.69097477,  -0.13459364,  0.71023828, -674.35431;
       0.011185364, -0.98438656, -0.17566407,  263.93491;
       0.72279227,  -0.11343517,  0.68169177, -677.69617;
       0,            0,           0,             1];
m1 = [-0.70978242,  -0.21733436, -0.67005581,  871.03137;
      -0.080154344, -0.92011875,  0.38334945, -266.29767;
      -0.69984591,   0.32580256,  0.63566369, -522.82904;
       0,            0,           0,             1];
m2 = [-0.97232783,   0.11579948, 0.20290148,  -28.091051;
      -0.020656202, -0.90772098, 0.41906551, -265.4436;
       0.23270552,   0.40327793, 0.88499433, -758.5166;
       0,            0,          0,             1];

% camera positions
c0 = m0(1:3,4);
c1 = m1(1:3,4);
c2 = m2(1:3,4);

ms = {m0, m1, m2};

% it has to be maximized
nW = arrayfun(@(nx,ny,nz) calc_measure(ms, [nx;ny;nz]), nx, ny, nz);

figure
hold on
slice(nx,ny,nz,nW,gx,gy,gz)
shading flat
scatter3(c0(1), c0(2), c0(3), 'r', 'filled');
scatter3(c1(1), c1(2), c1(3), 'g', 'filled');
scatter3(c2(1), c2(2), c2(3), 'b', 'filled');
%axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
hold off



function f = calc_measure(ms, p)
Ciws = cellfun(@(m) calc_Ci(m,p), ms, 'UniformOutput', false);
Ciw = sum(cat(3,Ciws{:}),3);
%Cw = inv(Ciw);
f = min(eig(Ciw));



% calculate the inverse of the covariance matrix
% for one camera and one observed point
% the result is given in the world coordinate system
function Ciw = calc_Ci(m, p)
R = m(1:3,1:3)';
t = -R*m(1:3,4);

Rt = [ R t ; zeros(1,3) 1];
pw = p;
pwh = [pw ; 1];
pch = Rt*pwh;
pc = pch(1:3);

% from ps3eye_intrinsics_red.xml (Avg_Reprojection_Error, Camera_Matrix)
e = 0.7758;
fx = 789.1510;
fy = 789.1510;
cx = 319.5;
cy = 239.5;

[Rotc x y d] = Rot3Dz2vect(pc); % rotation in the camera coord. system
global useDetectAngle
if useDetectAngle
    sigx = d * e / fx * cos(y)^2;
    sigy = d * e / fy * cos(x)^2;
else
    sigx = d * e / fx;
    sigy = d * e / fy;
end
six = sigx^(-2);
siy = sigy^(-2);
Ci = diag([six,siy,0]);

Rot = R' * Rotc;
Ciw = Rot*Ci*Rot';
