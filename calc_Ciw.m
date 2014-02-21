function [Ciw Ciws] = calc_Ciw(cams, pw)

% calculate the inverse of the resulting covariance matrix
% for one camera and one observed point
% the result is given in the world coordinate system
% cams: cell containing camera structs
% pw: point coordinate vector in the world coordinate system
% Ciw: inverse of the resulting covariance matrix
% Ciws: cell containing the inverses of the covariance matrices

% if only one camera is given, a cell is formulated
if isstruct(cams)
    cams = {cams};
end

Ciws = cellfun(@(cam) calc_Ciw_single(cam,pw), cams, 'UniformOutput', false);
Ciw = sum(cat(3,Ciws{:}),3);


function Ciw = calc_Ciw_single(cam, pw)

w2c = cam.w2c;
e = cam.e;
fx = cam.fx;
fy = cam.fy;
R = cam.R;

pwh = [pw ; 1];
pch = w2c * pwh;
pc = pch(1:end-1);

[Rotc x y d] = Rot3Dz2vect(pc); % rotation in the camera coord. system
global useDetectAngle
if useDetectAngle
    sigx = d * e / fx * cos(y)^2;
    sigy = d * e / fy * cos(x)^2;
else
    sigx = d * e / fx;
    sigy = d * e / fy;
end
six = sigx^(-2);
siy = sigy^(-2);
Ci = diag([six,siy,0]);

Rot = R' * Rotc;
Ciw = Rot*Ci*Rot';