function [] = main12()

% two camera observing the oobject at objX
% a third camera is placed
% calc_opt_plane applied

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


cams = {CreateCamera('pos', [   500;  1500 ]), ...
        CreateCamera('pos', [ -1000;   500 ]), ...
        CreateCamera('pos', [ -1000;  -500 ])};

polygon{1} = [  -500,  1500 ; ...
                -500,  2000 ; ...
                -300,  2000 ; ...
                 600,  1500 ; ...
                   0,  1100 ];
polygon{2} = [  -400,  1200 ; ...
                -700,  1600 ; ...
               -1100,   500 ; ...
                -300,  -300 ];
polygon{3} = [ -1200,   200 ; ...
                -700,  -700 ; ...
               -1000,  -600 ];
polygon = cellfun(@(p) [p;p(1,:)], polygon, 'uni', false);
num = numel(cams);
out = arrayfun(@(i) inpolygon(cams{i}.pos(1), cams{i}.pos(2), polygon{i}(:,1),polygon{i}(:,2)), 1:num);
if ~all(out)
    error(['camera number i is not placed in polygon with number i (i = ' num2str(find(~out)) ')'])
end

objX = [0;0];

            
hold on
axis('equal')
DrawCamera(cams)

cellfun(@(p) plot(p(:,1),p(:,2)), polygon); % plot all polygon
plot(objX(1), objX(2), 'g*') % plot objX

for i=1:6
    j = mod(i,num);
    if(j==0)
        j=3;
    end
    cams_i = cams;
    cams_i(j) = [];
    xopt = calc_opt_polygon(cams_i, polygon{j}, objX);
    cams{j} = CreateCamera('pos', xopt);
    plot(xopt(1), xopt(2), 'r*') % plot all xopts with red
end

hold off


