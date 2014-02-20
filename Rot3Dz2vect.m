function [R x y d] = Rot3Dz2vect(v)

% v: the z axis is rotated into the direction of v
%    rotation first around the x axis (elevation), second around the y axis
% x: rotation angle around the x axis
% y: rotation angle around the y axis
% d: distance, the length of the v vector
% R: 3x3 matrix
% usage: R * C * R'
% usage: R * p

% for checking
% rotate the coordinate system so
% that the image of the z become x,
% then rotate x to become v
% R = Rot3Dx2vect(v) * Rot3D('y',pi/2) * Rot3D('z',pi/2);

% misced coordinates
% new - old
%  x  -  x
%  y  -  -z
%  z  -  y
[y x d] = cart2sph(v(1), v(3), -v(2));
y = -(y-pi/2); % angle in the xz (originally xy) plane, but from the x axis
R = Rot3D('y',y) * Rot3D('x',x);
