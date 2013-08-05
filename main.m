
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



%%--- new ---
%the location of the new camera
[nX,nY] = meshgrid(90:1:130, 60:1:60);

ns1 = size(nX, 1);
ns2 = size(nX, 2);

nW = zeros(gs1,gs2,ns1,ns2);
ndW = zeros(gs1,gs2,ns1,ns2);




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


    %progress
    i,j
    fflush(stdout);

    %for every location of the new camera
    for m=1:ns1
      for n=1:ns2
        x = nX(m,n);
        y = nY(m,n);
        alpha = GetAlpha2D(dmX-x, dmY-y);
        saCov_n = CalculateCovariance([cam,CreateCamera(alpha,[x;y])], X);
        sCovRes_n = CombineGaussians(saCov_n);
        valid_n = sCovRes_n.valid;
        C_n = sCovRes_n.C;
        Ci_n = sCovRes_n.Ci;
        nW(i,j,m,n) = det(Ci_n);
        ndW(i,j,m,n) = nW(i,j,m,n) * dYX(i,j);
      end
    end

  end
end

hold off


%calculate the best location and orientation for the new camera
%best means, where it improves the most
MM1 = sum(ndW, 1);
MM2 = sum(MM1, 2);
[MM3,I3] = max(MM2, [], 3);
[MM4,I4] = max(MM3, [], 4);

m2=I4(1, 1, 1, 1);
m1=I3(1, 1, 1,m2);

mX = nX(m1,m2);
mY = nY(m1,m2);
mAlpha = GetAlpha2D(dmX-mX, dmY-mY);

figure(1)
hold on
DrawCamera(CreateCamera(mAlpha,[mX;mY]), "r");
hold off


%results
%max indexes
m1,m2
%location and orientation values
mX,mY,mAlpha


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

