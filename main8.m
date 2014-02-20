function [] = main8()

% calculate the measure in a volume
% using virtual cameras with parameters defined in ICSSE_2013 publication
% the measure is the smallest eigenvalue of the inverse cov. matr.
% it has to be maximized
% the figure is contourslice

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


gx = -1000:100:1000;
gy = 100:100:1000;
gz = -1000:200:1000;

[nx,ny,nz] = meshgrid(gx,gy,gz);

c0 = [-1000; 0; 0];
c1 = [1000; 0; 0];

m0 = [eye(3) c0 ; zeros(1,3) 1];
m1 = [eye(3) c1 ; zeros(1,3) 1];
ms = {m0, m1};


% it has to be maximized
nW = arrayfun(@(nx,ny,nz) calc_measure(ms, [nx;ny;nz]), nx, ny, nz);

figure
hold on
contourslice(nx,ny,nz,nW,[],[],gz)
shading flat
scatter3(c0(1), c0(2), c0(3), 'r', 'filled');
scatter3(c1(1), c1(2), c1(3), 'g', 'filled');
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
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
