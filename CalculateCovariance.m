function out = CalculateCovariance(cam, X)

%cam: camera struct
%X: 2x1 vector
%out: struct containing the fields: Ci, mu


out.Ci = zeros(2,2);
out.mu = X;

%field of view contains the point
for [val, key] = cam.normal_vectors
  if (X-cam.T)' * val <= 0
    return;
  end
end


%camera and point distance
t = norm(X-cam.T);

RT = [ cam.R  cam.T ;
       0  0     1   ] ;

X_w_h = [ X ;
          1 ] ;

X_c_h = RT * X_w_h;



alfa = atan( X_c_h(2) / X_c_h(1) );

sig = t * cos(alfa)^2 / cam.f_mm * cam.e_mm;

Ci = [ 0     0    ;
       0  1/sig^2 ] ;

%world <= camera <= Ci
Rot = cam.R * Rot2D(alfa)';

out.Ci = Rot * Ci * Rot';
out.mu = X;

