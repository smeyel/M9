
%file for trying the actual problem

clear
clc

%warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath

global useFoV
useFov=false;

X=95;
Y=0;

%%--- camera ---
%the location and oreientation of the cameras
x1=10;
y1=50;
cam(1) = CreateCamera(GetAlpha2D(X-x1, Y-y1), [x1;y1]);
x2=10;
y2=-50;
cam(2) = CreateCamera(GetAlpha2D(X-x2, Y-y2), [x2;y2]);


%%2D
figure(1); clf;
hold on
axis([0 150 -60 70], 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

DrawCamera(cam, 'b')

csCovRes = CalculateCovariance(cam, [X;Y]);
gsCovRes = CalculateResultingCovariance(cam, [X;Y]);
h = my_2D_error_ellipse(100*inv(csCovRes(1).Ci_kamu), [X;Y], 'conf', 0.95);
h = my_2D_error_ellipse(100*inv(csCovRes(2).Ci_kamu), [X;Y], 'conf', 0.95);
h = my_2D_error_ellipse(100*gsCovRes.C, [X;Y], 'conf', 0.95);
%plot([cam(1).pos(1) X], [cam(1).pos(2) Y])
%plot([cam(2).pos(1) X], [cam(2).pos(2) Y])

hold off




%save
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), 'figures/covariance_ellipses_one.eps')

