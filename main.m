
%file for trying the actual problem

clear
clc

%warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath

global useFoV;
useFoV=false;


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

rx = (max(max(gX)) - min(min(gX))) / size(gX,2) / 2;
ry = (max(max(gY)) - min(min(gY))) / size(gY,1) / 2;
dArea = [min(min(gX(dYX>0)))-rx min(min(gY(dYX>0)))-ry ; ...
         max(max(gX(dYX>0)))+rx min(min(gY(dYX>0)))-ry ; ...
         max(max(gX(dYX>0)))+rx max(max(gY(dYX>0)))+ry ; ...
         min(min(gX(dYX>0)))-rx max(max(gY(dYX>0)))+ry];

if sum(sum(dYX)) ~= 1
  error('dYX is not a valid density function!')
end

dmX = sum(sum(dYX .* gX));
dmY = sum(sum(dYX .* gY));




%%--- new ---
%the location of the new camera
[nX,nY] = meshgrid(85:1:135, 60:1:60);

nArea = [min(min(nX)) min(min(nY)) ; ...
         max(max(nX)) min(min(nY)) ; ...
         max(max(nX)) max(max(nY)) ; ...
         min(min(nX)) max(max(nY))];




%%--- wrong ---
wX = max(max(nX));
wY = min(min(nY));
wAlpha = GetAlpha2D(dmX-wX, dmY-wY);
wCam = CreateCamera(wAlpha,[wX;wY]);



%--- calculations ---
gsCovRes = arrayfun(@(gx,gy) CalculateResultingCovariance(cam, [gx;gy]), gX, gY);
gW = arrayfun(@(covres) det(covres.Ci), gsCovRes);

wsCovRes = arrayfun(@(gx,gy) CalculateResultingCovariance([cam,wCam], [gx;gy]), gX, gY);
wW = arrayfun(@(covres) det(covres.Ci), wsCovRes);
wdW = wW .* dYX;

ngX = repmat(reshape(gX, [size(gX), 1, 1]), [1, 1, size(nX)]);
ngY = repmat(reshape(gY, [size(gY), 1, 1]), [1, 1, size(nY)]);
ndYX = repmat(reshape(dYX, [size(dYX), 1, 1]), [1, 1, size(nX)]);
nnX = repmat(reshape(nX, [1, 1, size(nX)]), [size(gX), 1, 1]);
nnY = repmat(reshape(nY, [1, 1, size(nY)]), [size(gY), 1, 1]);
nsCovRes = arrayfun(@(gx,gy,nx,ny) CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(dmX-nx, dmY-ny),[nx;ny])], [gx;gy]), ...
    ngX, ngY, nnX, nnY);
nW = arrayfun(@(covres) det(covres.Ci), nsCovRes);
ndW = nW .* ndYX;
nsW = squeeze(sum(sum(ndW, 1), 2));



%calculate the best location and orientation for the new camera
%best means, where it improves the most
[MM1,I1] = max(nsW, [], 1);
[MM2,I2] = max(MM1, [], 2);
m2=I2(1, 1);
m1=I1(1, m2);

mX = nX(m1,m2);
mY = nY(m1,m2);
mAlpha = GetAlpha2D(dmX-mX, dmY-mY);
mCam = CreateCamera(mAlpha,[mX;mY]);


msCovRes = nsCovRes(:,:,m1,m2);
mW = nW(:,:,m1,m2);
mdW = ndW(:,:,m1,m2);


%results
%max indexes
m1,m2
%location and orientation values
mX,mY,mAlpha




%%2D
figure(1); clf;
hold on
axis([0 150 -60 70], 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
gv = [gsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    gsCovRes(gv), gX(gv), gY(gv));

hold off


%%2D
figure(2); clf;
hold on
axis([0 150 -60 70], 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
DrawCamera(mCam, 'g');

drawPolygon(dArea)
drawPolygon(nArea)

mv = [gsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    msCovRes(mv), gX(mv), gY(mv));

hold off




%%2D
figure(3); clf;
hold on
axis([0 150 -60 70], 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
DrawCamera(wCam, 'r');

drawPolygon(dArea)
drawPolygon(nArea)

wv = [gsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    wsCovRes(wv), gX(wv), gY(wv));

hold off




%%3D
figure(4); clf;
hold on
contour(gX,gY,gW, 60);
axis([0 150 -60 60], 'equal');
arrayfun(@(c) DrawCamera(c), cam);
hold off

figure(5); clf;
meshz(gX,gY,gW);
axis([0 150 -60 60]);
xlabel('x')
ylabel('y')


%save
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), 'figures/covariance_ellipses.eps')
saveas(figure(2), 'figures/covariance_ellipses_new.eps')
saveas(figure(3), 'figures/covariance_ellipses_wrong.eps')
saveas(figure(4), 'figures/contour.eps')
saveas(figure(5), 'figures/meshz.eps')

