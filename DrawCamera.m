function DrawCamera(cams, color)

% draw cameras
% cams: cell containing camera structs
% color: color string

if(nargin < 2)
    color='b';
end

% if only one camera is given, a cell is formulated
if isstruct(cams)
    cams = {cams};
end

cellfun(@(cam) DrawCamera_single(cam, color), cams);


function DrawCamera_single(cam, color)

global useDetectAngle;
l = 50; %length
c = cam.pos;
o = cam.ori;

if cam.dim == 2
    scatter(c(1), c(2), 'filled', color);
    if useDetectAngle
        plot([c(1) c(1)+l*o(1)], ...
             [c(2) c(2)+l*o(2)], color);
    end
else
    scatter3(c(1), c(2), c(3), 'filled', color);
    if useDetectAngle
        plot3([c(1) c(1)+l*o(1)], ...
              [c(2) c(2)+l*o(2)], ...
              [c(3) c(3)+l*o(3)], color);
    end
end


