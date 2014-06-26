function cam = CreateCamera(ori, pos)

%R: Rotation matrix, world <= camera, 2x2 matrix
%T: Translation matrix, world <= camera, 2x1 vector
%cam: camera struct


cam.ori = ori;
cam.pos = pos;

%%camera parameters
%camera <= world
cam.R = Rot2D(-ori);
cam.T = -cam.R * pos;
cam.RT = [ cam.R  cam.T ;
            0 0     1   ] ;


%cam.x_pixel_per_mm = 1;
%cam.y_pixel_per_mm = 1;
cam.pixel_per_mm = 1;

%cam.x_pixel = 640;
%cam.y_pixel = 480;
cam.pixel = 480;

%cam.x_f_pixel = 10;
%cam.y_f_pixel = 10;
cam.f_pixel = 400;

%cam.x_f_mm = cam.x_f_pixel / cam.x_pixel_per_mm;
%cam.y_f_mm = cam.y_f_pixel / cam.y_pixel_per_mm;
cam.f_mm = cam.f_pixel / cam.pixel_per_mm;

%cam.x_fov = atan(cam.x_pixel / cam.x_f_pixel);
%cam.y_fov = atan(cam.y_pixel / cam.y_f_pixel);
cam.fov = 2*atan((cam.pixel/2) / cam.f_pixel);

%cam.x_e_pixel = 1;
%cam.y_e_pixel = 1;
cam.e_pixel = 1;

%cam.x_e_mm = cam.x_e_pixel / cam.x_pixel_per_mm;
%cam.y_e_mm = cam.y_e_pixel / cam.y_pixel_per_mm;
cam.e_mm = cam.e_pixel / cam.pixel_per_mm;


%world <= camera
cam.normal_vectors.low = cam.R' * Rot2D(-cam.fov/2) * [0;1] ;
cam.normal_vectors.high = cam.R' * Rot2D(cam.fov/2) * [0;-1] ;
cam.normal_vectors.front = cam.R' * [1;0] ;

