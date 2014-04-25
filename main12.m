function [] = main12()

% two camera observing the oobject at objX
% a third camera is placed
% calc_opt_plane applied
% iterative placement of the 3 cameras, the optimum is not found

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


cams = {CreateCamera('pos', [  900;    50 ]), ...
        CreateCamera('pos', [ -1000;   400 ]), ...
        CreateCamera('pos', [  -900;  -500 ])};

polygon{1} = [   500,     0 ; ...
                700,   500 ; ...
                1100,  -100 ; ...
                 600,  -100 ; ...
                 600,  -500 ];
polygon{2} = [  0, 500/cos(pi/6) ; ...
                -500/cos(pi/3)+200, 0+200*tan(pi/6) ; ...
               -1100,   500 ; ...
                -800,   300 ];
polygon{3} = [  -500/cos(pi/3), 0 ; ...
                500*cos(-2*pi/3),  500*sin(-2*pi/3) ; ...
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


for i=1:15
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


camo = {CreateCamera('pos', [     500;   0 ]), ...
        CreateCamera('pos', [500*cos(2*pi/3); 500*sin(2*pi/3)]), ...
        CreateCamera('pos', [500*cos(-2*pi/3); 500*sin(-2*pi/3)])};
DrawCamera(camo, 'm')
Ci = calc_Ciw(cams, objX);
C = pinv(Ci);
min(eig(Ci))
Cio = calc_Ciw(camo, objX);
Co = pinv(Cio);
min(eig(Cio))

cd 'old'
my_2D_error_ellipse(1000*C,  objX, 'conf', 0.95, 'style', 'r');
my_2D_error_ellipse(1000*Co, objX, 'conf', 0.95);
cd '..'

hold off
