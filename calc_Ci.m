function Ciw = calc_Ci(cam, pw)

% calculate the inverse of the covariance matrix
% for one camera and one observed point
% the result is given in the world coordinate system
% cam: camera struct
% pw: point coordinate vector in the world coordinate system

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