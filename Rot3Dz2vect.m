function R = Rot3Dz2vect(v)

% v: the z axis is rotated into the direction of v
% R: 3x3 matrix
% usage: R * C * R'
% usage: R * p

if false
    R = Rot3Dx2vect(v) * Rot3D('y',pi/2) * Rot3D('z',pi/2);
else
    [y x d] = cart2sph(v(1), v(3), -v(2)); % misced coordinates
    y = -(y-pi/2); % angle in the xz (originally xy) plane, but from the x axis
    R = Rot3D('y',y) * Rot3D('x',x);
end
