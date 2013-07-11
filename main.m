
%file for trying the actual problem

%draw 2 cameras and the common covariance ellipse

clear
clc

warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath


%%--- camera ---
camsetup = 2;

if camsetup==1
  cam(1) = CreateCamera(0, [10;0]);
  cam(2) = CreateCamera(pi, [80;0]);
elseif camsetup==2
  cam(1) = CreateCamera(-pi/4, [10;50]);
  cam(2) = CreateCamera(pi/4, [10;-50]);
end



%%--- grid ---
gridsetup = 2;

if gridsetup==1
elseif gridsetup==2
  [gX,gY] = meshgrid(10:5:80, -20:2:20);
end
gr = size(gX, 1);
gc = size(gX, 2);

gZ = zeros(gr, gc);



%%2D
figure(1); clf;
hold on;
axis([-10 200 -100 100], "equal");

DrawCamera(cam)

for i=1:gr
  for j=1:gc
    X = [gX(i,j);gY(i,j)];
    saCov = CalculateCovariance(cam, X);
    sCovRes = CombineGaussians(saCov);
    C = sCovRes.C;
    Ci = sCovRes.Ci;
    if all(eig(Ci) > 1e-10) && all(eig(C) > 1e-10)
      h = my_2D_error_ellipse(C, X, 'conf', 0.95);
      gZ(i,j) = det(Ci);
    else
      gZ(i) = 0;
    end
  end
end

hold off


%%3D
figure(2); clf;
%plot3(gX(:), gY(:), gZ(:), ".");
surf(gX,gY,gZ);

