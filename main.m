
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
%gX(i,j) is the X with index j
%gY(i,j) is the Y with index i
[gX,gY] = meshgrid(20:10:140, -40:10:40);




%%--- density ---
%the density function
dYX = zeros(size(gX));

dYX(7:9,8:12) = 1/(3*5);

dArea = [gX(7,8)-5 gY(7,8)-5 ; gX(7,12)+5 gY(7,12)-5 ; gX(9,12)+5 gY(9,12)+5 ; gX(9,8)-5 gY(9,8)+5];

if sum(sum(dYX)) != 1
  error('dYX is not a valid density function!')
end

dmX = sum(sum(dYX .* gX));
dmY = sum(sum(dYX .* gY));




%%--- new ---
%the location of the new camera
[nX,nY] = meshgrid(90:1:130, 60:1:60);




%%--- wrong ---
wX = 60;
wY = 60;
wAlpha = GetAlpha2D(dmX-wX, dmY-wY);
wCam = CreateCamera(wAlpha,[wX;wY]);




%for every point in grid
for i=1:size(gX,1)
  for j=1:size(gX,2)

    X = gX(i,j);
    Y = gY(i,j);

    %resulting covariances
    gsCovRes(i,j) = CalculateResultingCovariance(cam, [X;Y]);
    gW(i,j) = det(gsCovRes(i,j).Ci);

    wsCovRes(i,j) = CalculateResultingCovariance([cam,wCam], [X;Y]);


    %progress
    i,j
    fflush(stdout);

    %for every location of the new camera
    for m=1:size(nX,1)
      for n=1:size(nX,2)
        x = nX(m,n);
        y = nY(m,n);
        alpha = GetAlpha2D(dmX-x, dmY-y);
        nsCovRes(i,j,m,n) = CalculateResultingCovariance([cam,CreateCamera(alpha,[x;y])], [X;Y]);
        nW(i,j,m,n) = det(nsCovRes(i,j,m,n).Ci);
        ndW(i,j,m,n) = nW(i,j,m,n) * dYX(i,j);
      end
    end

  end
end




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
mCam = CreateCamera(mAlpha,[mX;mY]);


msCovRes = nsCovRes(:,:,m1,m2);


%results
%max indexes
m1,m2
%location and orientation values
mX,mY,mAlpha




%%2D
figure(1); clf;
hold on
axis([0 150 -60 70], "equal");
xlabel("x")
ylabel("y", 'rotation', 0)

DrawCamera(cam)

%for every point in grid
for i=1:size(gX,1)
  for j=1:size(gX,2)
    if gsCovRes(i,j).valid
      h = my_2D_error_ellipse(10*gsCovRes(i,j).C, [gX(i,j);gY(i,j)], 'conf', 0.95);
    end
  end
end

hold off



%%2D
figure(2); clf;
hold on
axis([0 150 -60 70], "equal");
xlabel("x")
ylabel("y", 'rotation', 0)

DrawCamera(cam)
DrawCamera(mCam, "r");

drawPolygon(dArea)

%for every point in grid
for i=1:size(gX,1)
  for j=1:size(gX,2)
    if msCovRes(i,j).valid
      h = my_2D_error_ellipse(10*msCovRes(i,j).C, [gX(i,j);gY(i,j)], 'conf', 0.95);
    end
  end
end

hold off




%%2D
figure(3); clf;
hold on
axis([0 150 -60 70], "equal");
xlabel("x")
ylabel("y", 'rotation', 0)

DrawCamera(cam)
DrawCamera(wCam, "g");

drawPolygon(dArea)

%for every point in grid
for i=1:size(gX,1)
  for j=1:size(gX,2)
    if wsCovRes(i,j).valid
      h = my_2D_error_ellipse(10*wsCovRes(i,j).C, [gX(i,j);gY(i,j)], 'conf', 0.95);
    end
  end
end

hold off




%%3D
figure(4); clf;
hold on
contour(gX,gY,gW, 60);
axis([0 150 -60 60], "equal");
DrawCamera(cam)
hold off

figure(5); clf;
meshz(gX,gY,gW);
axis([0 150 -60 60]);
xlabel("x")
ylabel("y")


%save
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), "figures/covariance_ellipses.eps")
saveas(figure(2), "figures/covariance_ellipses_new.eps")
saveas(figure(3), "figures/covariance_ellipses_wrong.eps")
saveas(figure(4), "figures/contour.eps")
saveas(figure(5), "figures/meshz.eps")

