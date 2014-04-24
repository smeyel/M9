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

cams = {CreateCamera('pos', [-1000; 500; 0]), ...
        CreateCamera('pos', [-1000; -500; 0])};
    
polygon{1} = [ -500, 1500, 100 ; ...
                   -500, 2000, 100 ; ...
                    500, 1600, 100 ; ...
                    0, 1000, 100 ];
% polygon{2} = [ -400, 400 ; ...
%                    -400, 800 ; ...
%                     400, 1000 ; ...
%                     400, 400 ];
% polygon{3} = [ -600, 1500 ; ...
%                    1000, 500 ; ...
%                    -400, 600 ];
% polygon{4} = [ 1600, -2500 ; ...
%                    -2000, -1500 ; ...
%                    1400, -1800 ];
polygon = cellfun(@(p) [p;p(1,:)], polygon, 'uni', false);

objX = [0;0;0];

            
hold on
DrawCamera(cams)
% plot the first polygon and all of its vertices
%plot(polygon{1}(:,1),polygon{1}(:,2),...
%     polygon{1}(:,1),polygon{1}(:,2),'g*')

% plot all polygon and all vertices of all polygon
%cellfun(@(p) plot(p(:,1),p(:,2),...
%                  p(:,1),p(:,2),'g*'), polygon);

%polys = polygon;
%polys = cellfun(@(p) add_intersections(p,1), polys, 'uni', false);
%polys = cellfun(@(p) add_intersections(p,2), polys, 'uni', false);
% plot all polys and all vertices of all polys
%cellfun(@(p) plot(p(:,1),p(:,2),...
%                  p(:,1),p(:,2),'g*'), polys);

%plot(polygon{1}(:,1),polygon{1}(:,2)) % plot the first polygon
%xopt = calc_opt_polygon(cams, polygon{1}, objX);
%plot(xopt(1), xopt(2), 'r*') % plot xopt with red

%cellfun(@(p) plot(p(:,1),p(:,2)), polygon); % plot all polygon
cellfun(@(p) plot3(p(:,1),p(:,2),p(:,3)), polygon); % plot all polygon
xopts = cellfun(@(p) calc_opt_polygon(cams, p, objX), polygon, 'uni', false);
xopts = cell2mat(xopts)';
%plot(xopts(:,1), xopts(:,2), 'r*') % plot all xopts with red
plot3(xopts(:,1), xopts(:,2), xopts(:,3), 'r*') % plot all xopts with red

%plot(objX(1), objX(2), 'g*') % plot objX
plot3(objX(1), objX(2), objX(3), 'g*') % plot objX

hold off




return


function W = myfun(cams, x)
cams{length(cams)+1} = CreateCamera('pos', x);
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


function xopt = calc_opt_polygon(cams, polygon, objX)
% polygon: array of vertices

% Rotation and translation
[R ee] = eig(calc_Ciw(cams, objX));
T = objX;

% Transformation of the polygon to:
% object in origin and eigenvalues are parallel to the axes
polygon_cell = num2cell(polygon,2);
p_cell = cellfun(@(p) (R'*(p'-T))', polygon_cell, 'uni', false);
p = cell2mat(p_cell);

% Transformation of the cameras
cams_ori = cellfun(@(c) CreateCamera('pos', R'*(c.pos-T)), cams, 'uni', false);

% calc xopt and transform back
xopt = calc_opt_polygon_ori(cams_ori, p);
xopt = R*xopt+T;


function xopt = calc_opt_polygon_ori(cams, polygon)

if cams{1}.dim == 2
    ee = eig(calc_Ciw(cams, zeros(cams{1}.dim,1)));
    F = (cams{1}.e/cams{1}.fx)^(-2);
    E = abs(ee(2)-ee(1));
    D = ee(2)+ee(1);
    dist = sqrt(F/E);
    if ee(2)>ee(1)
        xmax = [0;dist];
    else
        xmax = [dist;0];
    end

    % xmax in polygon
    if inpolygon(xmax(1), xmax(2), polygon(:,1),polygon(:,2))
        xopt = xmax;
    % -xmax in polygon
    elseif inpolygon(-xmax(1), -xmax(2), polygon(:,1),polygon(:,2))
        xopt = -xmax;
    else
        % -xmax..xmax intersection with the polygon
        [XI,YI] = polyxpoly(polygon(:,1),polygon(:,2),[-xmax(1);xmax(1)],[-xmax(2);xmax(2)]);
        % intersection exists
        if ~isempty(XI)
            % the first intersection is selected
            xopt = [XI';YI'];
            xopt = xopt(:,1);
        % intersection does not exist
        else
            p2 = polygon;
            p2 = add_intersections(p2, 1);
            p2 = add_intersections(p2, 2);
            s2 = [p2(1:end-1,:) p2(2:end,:)]; % segments, each row: [x0 y0 x1 y1]

            s2_cell = num2cell(s2,2);
            [xopts qs] = cellfun(@(s) calc_opt_line_ori(cams, s), s2_cell, 'uni', false);
            [qopt qi] = min([qs{:}]);
            xopt = xopts{qi};
        end
    end

else

    % select three, non-collinear points
    p1 = polygon(1,:)';
    p2 = polygon(2,:)';
    for i=3:size(polygon,1)
        p = polygon(i,:)';
        if norm(cross(p-p1,p-p2)) ~= 0
            p3 = p;
            break
        end
    end
    if ~exist('p3')
        error('Points in polygon are collinear points!')
    end

    % generate base vectors
    w1 = p2-p1;
    w2 = p3-p1;
    v1 = w1 / norm(w1);
    w2n = w2 - v1'*w2*v1;
    v2 = w2n / norm(w2n);
    v3 = cross(v1,v2);

    % check coplanarity
    for i=1:size(polygon,1)
        if dot((polygon(i,:)'-p1),v3)
            error('Points in the polygon are not coplanar!')
        end
    end

    % find the global optimum (vopt) on the plane in the v1-v2 base
    vopt = fminunc(...
        @(v) -myfun(cams,p1+v(1)*v1+v(2)*v2), ... %fun
        [0;0], ... %x0
        []); %options
    % calculate xopt from vopt
    xopt = p1+vopt(1)*v1+vopt(2)*v2;

    % Transformation of the polygon and the optimum (xopt) to:
    % p1 in origin and v1,v2,v3 are parallel to the axes respectively
    R = [v1 v2 v3];
    T = p1;
    polygon_cell = num2cell(polygon,2);
    p_cell = cellfun(@(p) (R'*(p'-T))', polygon_cell, 'uni', false);
    p_plane = cell2mat(p_cell);
    x_plane = R'*(xopt-T);

    if inpolygon(x_plane(1), x_plane(2), p_plane(:,1),p_plane(:,2))
        return
    else
        p2 = polygon;
        p2 = add_intersections(p2, 1);
        p2 = add_intersections(p2, 2);
        p2 = add_intersections(p2, 3);
        s2 = [p2(1:end-1,:) p2(2:end,:)]; % segments, each row: [x0 y0 z0 x1 y1 z1]

        s2_cell = num2cell(s2,2);
        [xopts qs] = cellfun(@(s) calc_opt_line_ori(cams, s), s2_cell, 'uni', false);
        [qopt qi] = min([qs{:}]);
        xopt = xopts{qi};
    end

end


function [xopt q] = calc_opt_line_ori(cams, segment)
if cams{1}.dim == 3
    x0 = segment(1);
    y0 = segment(2);
    z0 = segment(3);
    x1 = segment(4);
    y1 = segment(5);
    z1 = segment(6);
    dx = x1-x0;
    dy = y1-y0;
    dz = z1-z0;
    [topt q] = fmincon(...
            @(t) -myfun(cams,[x0+t*dx;y0+t*dy;z0+t*dz]), ... %fun
            0, ... %x0
            [], [], ... %A, b
            [], [], ... %Aeq, beq
            0, ... %lb
            1, ... %ub
            [], ... %nonlcon
            []); %options
    xopt = [x0+topt*dx;y0+topt*dy;z0+topt*dz];
else
    x0 = segment(1);
    y0 = segment(2);
    x1 = segment(3);
    y1 = segment(4);
    dx = x1-x0;
    dy = y1-y0;
    [topt q] = fmincon(...
            @(t) -myfun(cams,[x0+t*dx;y0+t*dy]), ... %fun
            0, ... %x0
            [], [], ... %A, b
            [], [], ... %Aeq, beq
            0, ... %lb
            1, ... %ub
            [], ... %nonlcon
            []); %options
    xopt = [x0+topt*dx;y0+topt*dy];
end


function p_added = add_intersections(polygon, index)
% add the intersections to the polygon
% valid intersection:
%   the coordinate at the index changes the sign between -1 and +1
% polygon: array of vertices

spa = size(polygon,1)-1;
p_added = []; % empty
% for every edge in the polygon
for n=1:spa
    p_added = [p_added ; polygon(n,:)]; % add begin vertex
    i0 = polygon(n,index);
    i1 = polygon(n+1,index);
    % if valid intersection, signs are: -1 and +1
    if(i0*i1 < 0)
        l0 = abs(i0);
        l1 = abs(i1);
        x0 = polygon(n,:);
        x1 = polygon(n+1,:);
        new = x0 * l1/(l0+l1) + x1 * l0/(l0+l1);
        p_added = [p_added ; new];
    end
end
p_added = [p_added ; polygon(spa+1,:)]; % add end vertex
