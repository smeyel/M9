function cam = CreateCamera(varargin)

% Create camera model
% modes:
%   CreateCamera('oripos', ori, pos)
%   CreateCamera('c2w', c2w)
%   CreateCamera('w2c', w2c)
% Some camera parameters from: ps3eye_intrinsics_red.xml

% w2c: world => camera
% w2c = inv(m) = [ R t ;
%                  0 1]
% c2w: camera => world
% c2w = inv(Rt) = [ R' -R'*t ;
%                   0    1   ]
% campos = -R'*t;
% camori = R' * [0;0;1];

cam = struct(...
  'R', 0, ...
  't', 0, ...
  'ori', 0, ...
  'pos', 0, ...
  'c2w', 0, ...
  'w2c', 0);

if length(varargin) < 1
    error('Missing parameter!!!')
end

ii = 1;
while ii <= length(varargin)
    if strcmp(varargin{ii}, 'oripos')
        ori = varargin{ii+1};
        pos = varargin{ii+2};
        ii = ii + 3;
        dim = length(pos);
        validate_dim(dim);
        if dim == 3
            R = Rot3Dz2vect(ori)';
        else
            R = Rot3Dy2vect([ori;0])';
            R = R(1:dim,1:dim);
        end
        t = -R * pos;
        w2c = [R t ; zeros(1,dim) 1];
        c2w = [R' pos ; zeros(1,dim) 1];
    elseif strcmp(varargin{ii}, 'w2c')
        w2c = varargin{ii+1};
        ii = ii + 2;
        dim = max(size(w2c))-1;
        validate_dim(dim);
        R = w2c(1:dim,1:dim);
        t = w2c(1:dim,dim+1);
        pos = -R' * t;
        ori = R' * [zeros(dim-1,1);1];
        c2w = [ R' pos ; zeros(1,dim) 1];
	elseif strcmp(varargin{ii}, 'c2w')
        c2w = varargin{ii+1};
        ii = ii + 2;
        dim = max(size(c2w))-1;
        validate_dim(dim);
        R = c2w(1:dim,1:dim)';
        pos = c2w(1:dim,dim+1);
        t = -R * pos;
        ori = R' * [zeros(dim-1,1);1];
        w2c = [R t ; zeros(1,dim) 1];
    else
        error('No such mode!!!')
    end
end

cam.ori = ori/norm(ori);
cam.pos = pos;
cam.R = R;
cam.t = t;
cam.c2w = c2w;
cam.w2c = w2c;

% get camera parameters from ps3eye_intrinsics_red.xml
data = getData('ps3eye_intrinsics_red.xml');
dataNames = fieldnames(data);
for i = 1:numel(dataNames) 
    cam.(dataNames{i}) = data.(dataNames{i});
end


% some matrices from the camera data for the distortion calculation
if dim == 3
    cam.A = [ cam.fx    0     cam.cx ; ...
                0     cam.fy  cam.cy ; ...
                0       0       1    ];
else
    cam.A = [ cam.fx  cam.cx ; ...
                0       1 ];
end
cam.iR = inv(cam.A * cam.R);
cam.ir = reshape(cam.iR', 1, numel(cam.iR));



function validate_dim(dim)
if dim < 2 || 3 < dim
    error('Invalid dimension count!!!')
end
