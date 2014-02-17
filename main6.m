function [] = main6()

% calculate the resulting covariance matrix
% in case of one specific observed point
% with the cameras and measurement setup
% in the ICSSE_2013 publication

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



% camera positions
% from the ICSSE_2013 publication
% from the log
t0 = [-674.3543 ; 263.9349 ; -677.6962];
t1 = [871.0314 ; -266.2977 ; -522.8290];
t2 = [-28.0911 ; -265.4436 ; -758.5166];

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
        Ciw0 = calc_Ci(t0, p);
        Ciw1 = calc_Ci(t1, p);
        Ciw2 = calc_Ci(t2, p);
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
function Ciw = calc_Ci(t, p)

e = 0.7758; % from ps3eye_intrinsics_red.xml (Avg_Reprojection_Error)
% from the ICSSE_2013 publication, from the C++ code
fx = 789.1510;
fy = 789.1510;

v = p-t;
x = 0; % fx = fy = f, symmetric, it can be zero
[z y d] = cart2sph(v(1), v(2), v(3));
y = -y; % elevation from xy plane, but clockwise around y axis

%  cam -> world
%      ->  x
%   x  ->  y
%   y  ->  z
sigy = d * e / fx;
sigz = d * e / fy;
siy = sigy^(-2);
siz = sigz^(-2);
Ci = diag([0,siy,siz]);


Rz = [cos(z) -sin(z) 0 ; sin(z) cos(z) 0 ; 0 0 1];
Ry = [cos(y) 0 sin(y) ; 0 1 0 ; -sin(y) 0 cos(y)];
Rx = [1 0 0 ; 0 cos(x) -sin(x) ; 0 sin(x) cos(x)];
Ciw = Rz*Ry*Rx*Ci*Rx'*Ry'*Rz';
