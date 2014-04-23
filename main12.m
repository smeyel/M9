function [] = main12()

% two camera observing the origin
% a third camera is placed
% fmincon in 2D
% calc_opt_plane: translation and rotation applied (on object, polygon)
%   to get the default case (origin and eigenvectors parallel to the axes)
%   2D OK, in the 3D case correction needed, TODO

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath

cams = {CreateCamera('oripos', [0;1], [-1000; 500]), ...
        CreateCamera('oripos', [0;1], [-1000; -500])};
    
polygon{1} = [ -500, 500 ; ...
                   -500, 2000 ; ...
                    500, 2000 ; ...
                    500, 500 ];
polygon{2} = [ -400, 400 ; ...
                   -400, 1000 ; ...
                    400, 1000 ; ...
                    400, 400 ];
polygon{3} = [ -600, 1500 ; ...
                   1000, 500 ; ...
                   -400, 600 ];
polygon = cellfun(@(p) [p;p(1,:)], polygon, 'uni', false);

objX = [0;0];

            
hold on
DrawCamera(cams)
%cellfun(@(p) plot(p(:,1),p(:,2)), polygon); % plot all polygon
plot(polygon{1}(:,1),polygon{1}(:,2)) % plot the first polygon
xopts = calc_opt_plane(cams, polygon{1}, objX);
plot(xopts(:,1), xopts(:,2), 'r*'); % plot all xopts with red
hold off




return


function W = myfun(cams, x)
cams{length(cams)+1} = CreateCamera('oripos', [0;1], x);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));



function stop = outputfun(x, optimValues, state)
% output function for the solver
% plots the points and connects the consecutive ones
global firstCall
global pre
if ~firstCall
    plot([x(1) pre(1)], [x(2) pre(2)], 'k')
end
firstCall = false;
plot(x(1), x(2), 'g*')
pre = x;
stop = 0;


function [xopts ins] = calc_opt_plane(cams, polygon, objX)
% polygon: array of vertices

% Rotation and translation
[R ee] = eig(calc_Ciw(cams, objX));
ee = diag(ee);
T = objX;

% Transformation of the polygon to:
% object in origin and eigenvalues are parallel to the axes
polygon_cell = num2cell(polygon,2);
p_cell = cellfun(@(p) (R'*(p'-T))', polygon_cell, 'uni', false);
p = cell2mat(p_cell);

% Cutting the polygon
% polys: cell of polygons
polys = {p};

ss = numel(polys);
dim = numel(objX);
xopts = zeros(ss,dim);
ins = zeros(ss);
for n = 1:ss
    p = polys{n};

    if cams{1}.dim == 2
        F = (cams{1}.e/cams{1}.fx)^(-2);
        E = abs(ee(2)-ee(1));
        D = ee(2)+ee(1);
        dist = sqrt(F/E);
        if ee(2)>ee(1)
            xmax = [0;dist];
        else
            xmax = [dist;0];
        end

        if inpolygon(xmax(1), xmax(2), p(:,1),p(:,2))
            xopt = xmax;
            in = true;
        else
            [XI,YI] = polyxpoly(p(:,1),p(:,2),[0;xmax(1)],[0;xmax(2)]);
            if isempty(XI)
                xopt = xmax;
                in = false;
            else
                % the first intersection is selected
                xopt = [XI';YI'];
                xopt = xopt(:,1);
                in = true;
            end
        end

        xopt = R*xopt+T;

    else

        % TODO: only moved, correction needed
        start = [500; 1500];
        lb = [0;10];
        xopt = fmincon(...
            @(v) -myfun(cams,v), ... %fun
            start, ... %x0
            [], [], ... %A, b
            [], [], ... %Aeq, beq
            lb, ... %lb
            [], ... %ub
            [], ... %nonlcon
            optimset('OutputFcn', @outputfun)); %options
        plot(xopt(1), xopt(2), 'r*')
        in = inpolygon(xopt(1), xopt(2), p(:,1),p(:,2));

    end

    xopts(n,:) = xopt';
    ins(n) = in;

end

