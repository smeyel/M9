function [R y z d] = Rot3Dx2vect(v)

% v: the x axis is rotated into the direction of v
%    rotation first around the y axis (elevation), second around the z axis
% y: rotation angle around the y axis
% z: rotation angle around the z axis
% d: distance, the length of the v vector
% R: 3x3 matrix
% usage: R * C * R'
% usage: R * p

[z y d] = cart2sph(v(1), v(2), v(3));
y = -y; % elevation from xy plane, but clockwise around y axis
R = Rot3D('z',z) * Rot3D('y',y);
