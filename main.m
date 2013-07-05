
%file for trying the actual problem

%draw 2 cameras and the common covariance ellipse

myAddPath

figure(1); clf;
hold on;
axis([-10 200 -100 100], "equal");


c1 = CreateCamera();
c2 = CreateCamera(Rot2D(pi/4), [10;-50]);
DrawCamera(c1)
DrawCamera(c2)

X = [75;6];

C1 = CalculateCovariance(c1, X);
C2 = CalculateCovariance(c2, X);
C = CombineGaussians([C1;C2]).C;
h = my_2D_error_ellipse(C, X, 'conf', 0.95);


hold off

