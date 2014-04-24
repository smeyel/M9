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

cams = {CreateCamera('pos', [-1000; 500; 0]), ...
        CreateCamera('pos', [-1000; -500; 0])};
    
polygon{1} = [ -500, 1500, 100 ; ...
                   -500, 2000, 100 ; ...
                    500, 1600, 100 ; ...
                    0, 1000, 100 ];
polygon = cellfun(@(p) [p;p(1,:)], polygon, 'uni', false);

objX = [0;0;0];

            
hold on
DrawCamera(cams)

%cellfun(@(p) plot(p(:,1),p(:,2)), polygon); % plot all polygon
cellfun(@(p) plot3(p(:,1),p(:,2),p(:,3)), polygon); % plot all polygon
xopts = cellfun(@(p) calc_opt_polygon(cams, p, objX), polygon, 'uni', false);
xopts = cell2mat(xopts)';
%plot(xopts(:,1), xopts(:,2), 'r*') % plot all xopts with red
plot3(xopts(:,1), xopts(:,2), xopts(:,3), 'r*') % plot all xopts with red

%plot(objX(1), objX(2), 'g*') % plot objX
plot3(objX(1), objX(2), objX(3), 'g*') % plot objX

hold off


