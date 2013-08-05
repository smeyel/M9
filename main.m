
%file for trying the actual problem

clear
clc

warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath

global useFoV=false;


%%--- camera ---
%the location and oreientation of the cameras
cam(1) = CreateCamera(-pi/4, [10;50]);
cam(2) = CreateCamera(pi/4, [10;-50]);



%%--- grid ---
%the accuracy is calculated in this points
[gX,gY] = meshgrid(20:10:140, -40:10:40);

gs1 = size(gX, 1);
gs2 = size(gX, 2);

gW = zeros(gs1, gs2);



%%--- density ---
ds1 = gs1;
ds2 = gs2;
dYX = zeros(ds1, ds2);

dYX(7:9,8:12) = 1/(3*5);

dmX = sum(sum(dYX .* gX));
dmY = sum(sum(dYX .* gY));

dW = zeros(ds1, ds2);




%%2D
figure(1); clf;
hold on;
axis([0 150 -60 60], "equal");
xlabel("x")
ylabel("y", 'rotation', 0)

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
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), "figures/covariance_ellipses.eps")
saveas(figure(2), "figures/contour.eps")
saveas(figure(3), "figures/meshz.eps")

