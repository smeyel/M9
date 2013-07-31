
%file for trying the actual problem

clear
clc

warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath

global useFoV=false;


%%--- camera ---
%the location and oreientation of the cameras
camsetup = 2;

if camsetup==1
  cam(1) = CreateCamera(0, [10;0]);
  cam(2) = CreateCamera(pi, [80;0]);
elseif camsetup==2
  cam(1) = CreateCamera(-pi/4, [10;50]);
  cam(2) = CreateCamera(pi/4, [10;-50]);
end



%%--- grid ---
%the accuracy is calculated in this points
gridsetup = 2;

if gridsetup==1
elseif gridsetup==2
  [gX,gY] = meshgrid(20:10:140, -40:10:40);
end
gs1 = size(gX, 1);
gs2 = size(gX, 2);

gW = zeros(gs1, gs2);




%%2D
figure(1); clf;
hold on;
axis([0 150 -60 60], "equal");

DrawCamera(cam, "k")

%for every point in grid
for i=1:gs1
  for j=1:gs2

    %covariance ellipses with cameras in the camsetup are drawn
    X = [gX(i,j);gY(i,j)];
    saCov = CalculateCovariance(cam, X);
    sCovRes = CombineGaussians(saCov);
    valid = sCovRes.valid;
    C = sCovRes.C;
    Ci = sCovRes.Ci;
    if valid
      h = my_2D_error_ellipse(10*C, X, 'conf', 0.95, 'style', "k");
    end
    gW(i,j) = det(Ci);
  end
end

hold off


%%3D
figure(2); clf;
hold on
colormap([0 0 0])
contour(gX,gY,gW, 60);
axis([0 150 -60 60], "equal");
DrawCamera(cam, "k")
hold off

figure(3); clf;
meshz(gX,gY,gW);
axis([0 150 -60 60]);
xlabel("x")
ylabel("y")


%save
saveas(figure(1), "figures/covariance_ellipses.eps")
saveas(figure(2), "figures/contour.eps")
saveas(figure(3), "figures/meshz.eps")

