function [] = main6()

% the standard deviation of marker location estimations,
% shown for the 3 axes separately
% the cameras and measurement setup is defined
% in the ICSSE_2013 publication

% Rt: world => camera
% Rt = inv(m) = [ R t ;
%                 0 1]
% m: camera => world
% m = inv(Rt) = [ R' -R'*t ;
%                 0    1   ]
% campos = -R'*t;
% camori = R' * [0;0;1];

close all
clear
clc

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

means = LocationMean3Ray;
cc = max(size(means)); % column count
TheoryLocationStd = zeros(cc,3);
for i=1:cc

    % observed point
    p = means(i,:)';

    if any(isnan(p))
        sigres = [0;0;0];
    else
        % covariance matrix inverses and the resulting covariace matrix
        Ciw0 = calc_Ci(m0, p);
        Ciw1 = calc_Ci(m1, p);
        Ciw2 = calc_Ci(m2, p);
        Ciw = Ciw0 + Ciw1 + Ciw2;
        Cw = inv(Ciw);

        % variance in the x-y-z directions
        [V D] = eig(Cw);
        sigres = V.^2 * diag(D);
    end
    TheoryLocationStd(i,1:3) = sigres';

end

bar(TheoryLocationStd(:,:));
axis([1 size(TheoryLocationStd,1) 0 5]);
legend('X','Y','Z');
xlabel('Location index');
ylabel('Standard deviation');

return


% figure
figure
hold on
error_ellipse(Cw, p, 'conf', 0.95);
scatter3(t0(1),t0(2),t0(3), 'r', 'filled');
scatter3(t1(1),t1(2),t1(3), 'g', 'filled');
scatter3(t2(1),t2(2),t2(3), 'b', 'filled');
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');
hold off


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

d = norm(pc);
sigx = d * e / fx;
sigy = d * e / fy;
six = sigx^(-2);
siy = sigy^(-2);
Ci = diag([six,siy,0]);

Rot = R' * Rot3Dz2vect(pc);
Ciw = Rot*Ci*Rot';

