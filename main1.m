function [] = main1()

% file for trying the actual problem
% Display the covariance ellipses for two cameras
% Display the resulting covariance ellipse

%% notes
% single letters in variable names:
% e - ellipse
% c - camera

%% preparation

close all
clear
clc

myAddPath

%colormap([zeros(63,3) ; ones(1,3)]);
%warning('off', 'Octave:possible-matlab-short-circuit-operator');

global useFoV;
useFoV=false;


%% local variable definitions

% --- displayArea ---
displayArea = [0 150 -60 70];


% --- ellipse center ---
eX=95;
eY=0;

% --- camera ---
%the location and orientation of the cameras
cX1=10;
cY1=50;
cX2=10;
cY2=-50;
cam(1) = CreateCamera(cart2pol(eX-cX1, eY-cY1), [cX1;cY1]);
cam(2) = CreateCamera(cart2pol(eX-cX2, eY-cY2), [cX2;cY2]);


%% Covariance ellipses

csCovRes = arrayfun(@(camera) CalculateCovariance(camera, [eX;eY]), cam);
esCovRes = CalculateResultingCovariance(cam, [eX;eY]);

%figure
fig_covariance_ellipses_one = figure; clf;
hold on
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);

my_2D_error_ellipse(100*pinv(csCovRes(1).Ci_kamu), [eX;eY], 'conf', 0.95);
my_2D_error_ellipse(100*pinv(csCovRes(2).Ci_kamu), [eX;eY], 'conf', 0.95);
my_2D_error_ellipse(100*esCovRes.C, [eX;eY], 'conf', 0.95);
%plot([cam(1).pos(1) eX], [cam(1).pos(2) eY])
%plot([cam(2).pos(1) eX], [cam(2).pos(2) eY])

hold off


%% save figures
saveas(fig_covariance_ellipses_one, 'figures/covariance_ellipses_one.eps')
