function [R x z d] = Rot3Dy2vect(v)

% v: the x axis is rotated into the direction of v
%    rotation first around the x axis (elevation), second around the z axis
% x: rotation angle around the x axis
% z: rotation angle around the z axis
% d: distance, the length of the v vector
% R: 3x3 matrix
% usage: R * C * R'
% usage: R * p

[z x d] = cart2sph(v(1), v(2), v(3));
z = z-pi/2; % angle in the xy plane, but from the x axis
R = Rot3D('z',z) * Rot3D('x',x);
