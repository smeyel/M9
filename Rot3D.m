function R = Rot3D(axis, angle)

% axis: 'x' or 'y' or 'z'
% angle
% R: 3x3 matrix
% usage: R * C * R'
% usage: R * p

if strcmpi(axis, 'x')
    R = [1 0 0 ; 0 cos(angle) -sin(angle) ; 0 sin(angle) cos(angle)];

elseif strcmpi(axis, 'y')
    R = [cos(angle) 0 sin(angle) ; 0 1 0 ; -sin(angle) 0 cos(angle)];

elseif strcmpi(axis, 'z')
    R = [cos(angle) -sin(angle) 0 ; sin(angle) cos(angle) 0 ; 0 0 1];

else
    error(['No axis named: ' axis])
end
