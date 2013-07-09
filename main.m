
%file for trying the actual problem

%draw 2 cameras and the common covariance ellipse

clear
clc

myAddPath

figure(1); clf;
hold on;
axis([-10 200 -100 100], "equal");


cam(1) = CreateCamera();
cam(2) = CreateCamera(pi/4, [10;-50]);
DrawCamera(cam)

X = [75;6];

Ci_mu = CalculateCovariance(cam, X);
C = CombineGaussians(Ci_mu).C;
h = my_2D_error_ellipse(C, X, 'conf', 0.95);


hold off

