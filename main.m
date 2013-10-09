function [] = main()

% file for trying the actual problem
% Wellnes with specific camera configuration
% determination of the position of the new camera
% optimal and not optimal placement of the new camera


%% notes

% single letters in variable names:
% g - grid
% s - struct
% d - density
% r - reshaped
% n - new
% m - max (optimal)
% w - wring (not optimal)
% v - valid


%% preparation

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


% --- camera ---
%the location and oreientation of the cameras
cam(1) = CreateCamera(-pi/4, [10;10]);
cam(2) = CreateCamera(pi/4, [10;-10]);


% --- grid ---
%the localization accuracy is calculated in this points
%gX(i,j) is the X with index j
%gY(i,j) is the Y with index i
[gX,gY] = meshgrid(40:10:60, -10:10:10);

%grid step
gsx = (max(max(gX)) - min(min(gX))) / (size(gX,2)-1);
gsy = (max(max(gY)) - min(min(gY))) / (size(gY,1)-1);


% --- density ---
%the density function
dYX = zeros(size(gX));

% uniform distribution
dYX(:) = 1/(size(gX,1)*size(gX,2));

if sum(sum(dYX)) ~= 1
  error('dYX is not a valid density function!')
end

%the expectation (mean) of the density function
dmX = sum(sum(dYX .* gX));
dmY = sum(sum(dYX .* gY));

%area of the density function, displayed
%rectangle increased with the step size
dArea = [min(min(gX(dYX>0)))-gsx/2 min(min(gY(dYX>0)))-gsy/2 ; ...
         max(max(gX(dYX>0)))+gsx/2 min(min(gY(dYX>0)))-gsy/2 ; ...
         max(max(gX(dYX>0)))+gsx/2 max(max(gY(dYX>0)))+gsy/2 ; ...
         min(min(gX(dYX>0)))-gsx/2 max(max(gY(dYX>0)))+gsy/2];


% --- new ---
%the location of the new camera
[nX,nY] = meshgrid(85:1:100, -50:1:50);

%area of the new camera placement, displayed
%rectangle
nArea = [min(min(nX)) min(min(nY)) ; ...
         max(max(nX)) min(min(nY)) ; ...
         max(max(nX)) max(max(nY)) ; ...
         min(min(nX)) max(max(nY))];


% --- wrong ---
%wrong camera placement, not optimal position
%at one of the corners of the new camera rectangle
wX = max(max(nX));
wY = min(min(nY));
wAlpha = GetAlpha2D(dmX-wX, dmY-wY);
wCam = CreateCamera(wAlpha,[wX;wY]);





%% Wellness on the grid

gsCovRes = arrayfun(@(gx,gy) CalculateResultingCovariance(cam, [gx;gy]), gX, gY);
gW = arrayfun(@(covres) det(covres.Ci), gsCovRes);

figure(1); clf;
hold on
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
gv = [gsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    gsCovRes(gv), gX(gv), gY(gv));

hold off


%% Wellness with the new camera

%reshaped for 4D (gx-gy-nx-ny)
rgX = repmat(reshape(gX, [size(gX), 1, 1]), [1, 1, size(nX)]);
rgY = repmat(reshape(gY, [size(gY), 1, 1]), [1, 1, size(nY)]);
rdYX = repmat(reshape(dYX, [size(dYX), 1, 1]), [1, 1, size(nX)]);
rnX = repmat(reshape(nX, [1, 1, size(nX)]), [size(gX), 1, 1]);
rnY = repmat(reshape(nY, [1, 1, size(nY)]), [size(gY), 1, 1]);

rngsCovRes = arrayfun(@(gx,gy,nx,ny) CalculateResultingCovariance([cam,CreateCamera(GetAlpha2D(dmX-nx, dmY-ny),[nx;ny])], [gx;gy]), ...
    rgX, rgY, rnX, rnY);
rngW = arrayfun(@(covres) det(covres.Ci), rngsCovRes);
rngdW = rngW .* rdYX;
nW = squeeze(sum(sum(rngdW, 1), 2));

figure(4); clf;
hold on
contour(nX,nY,nW, 60);
axis(displayArea, 'equal');
arrayfun(@(c) DrawCamera(c), cam);
hold off

figure(5); clf;
meshz(nX,nY,nW);
axis(displayArea);
xlabel('x')
ylabel('y')


%% Wellness with the max (optimally placed) camera

%calculate the best location and orientation for the new camera
%best means, where it improves the most
[MM1,I1] = max(nW, [], 1);
[MM2,I2] = max(MM1, [], 2);
m2=I2(1, 1);
m1=I1(1, m2);

mX = nX(m1,m2);
mY = nY(m1,m2);
mAlpha = GetAlpha2D(dmX-mX, dmY-mY);
mCam = CreateCamera(mAlpha,[mX;mY]);

mgsCovRes = rngsCovRes(:,:,m1,m2);
mgW = rngW(:,:,m1,m2);
mgdW = rngdW(:,:,m1,m2);

figure(2); clf;
hold on
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
DrawCamera(mCam, 'g');

drawPolygon(dArea)
drawPolygon(nArea)

mv = [mgsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    mgsCovRes(mv), gX(mv), gY(mv));

hold off


%% Wellness with the wrong (not optimally placed) camera

wgsCovRes = arrayfun(@(gx,gy) CalculateResultingCovariance([cam,wCam], [gx;gy]), gX, gY);
wgW = arrayfun(@(covres) det(covres.Ci), wgsCovRes);
wgdW = wgW .* dYX;

figure(3); clf;
hold on
axis(displayArea, 'equal');
xlabel('x')
ylabel('y', 'rotation', 0)

arrayfun(@(c) DrawCamera(c), cam);
DrawCamera(wCam, 'r');

drawPolygon(dArea)
drawPolygon(nArea)

wv = [wgsCovRes.valid];
arrayfun(@(covres, gx, gy) my_2D_error_ellipse(10*covres.C, [gx;gy], 'conf', 0.95), ...
    wgsCovRes(wv), gX(wv), gY(wv));

hold off


%% save figures
saveas(figure(1), 'figures/covariance_ellipses.eps')
saveas(figure(2), 'figures/covariance_ellipses_new.eps')
saveas(figure(3), 'figures/covariance_ellipses_wrong.eps')
saveas(figure(4), 'figures/contour.eps')
saveas(figure(5), 'figures/meshz.eps')
