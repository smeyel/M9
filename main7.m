function [] = main7()

% calculate the measure
% in the volume defined in the ICSSE_2013 publication
% the measure is the smallest eigenvalue of the inverse cov. matr.
% it has to be maximized
% the cameras and measurement setup is defined in ICSSE_2013 publication
% the figure is slice

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


data = getData('M1R1_stats');
LocationMeanAll = data.LocationMeanAll;


minLoc = min(LocationMeanAll)';
maxLoc = max(LocationMeanAll)';

step = 10;
gx = minLoc(1):step:maxLoc(1);
gy = minLoc(2):step:maxLoc(2);
gz = minLoc(3):step:maxLoc(3);
[nx,ny,nz] = meshgrid(gx,gy,gz);


data = getData('M1R1_c2w');
dataNames = fieldnames(data);
for i = 1:numel(dataNames)
    cams{i} = CreateCamera('c2w', data.(dataNames{i}));
end


% measure, it has to be maximized
nW = arrayfun(@(nx,ny,nz) min(eig(calc_Ciw(cams, [nx;ny;nz]))), nx, ny, nz);

figure
hold on
slice(nx,ny,nz,nW,gx,gy,gz)
shading flat
DrawCamera(cams{1}, 'r');
DrawCamera(cams{2}, 'g');
DrawCamera(cams{3}, 'b');
grid on
axis('equal')
view(-165,-65)
xlim([-700,1000])
ylim([-300,300])
zlim([-800,300])
xlabel('x');
ylabel('y');
zlabel('z');
hold off
