function R = Rot3D(ax, ay, az)

% ax, ay, az: rotation angles around the axes, respectively
% R: 3x3 matrix
% usage: R * C * R'

Rz = [cos(az) -sin(az) 0 ; sin(az) cos(az) 0 ; 0 0 1];
Ry = [cos(ay) 0 sin(ay) ; 0 1 0 ; -sin(ay) 0 cos(ay)];
Rx = [1 0 0 ; 0 cos(ax) -sin(ax) ; 0 sin(ax) cos(ax)];

R = Rz*Ry*Rx;
