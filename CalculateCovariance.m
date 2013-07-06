function out = CalculateCovariance(cam, X)

%cam: nx1 struct array of camera structs
%X: 2x1 vector
%out: nx1 struct array containing the fields: Ci, mu


for i = 1:length(cam)

  out(i).Ci = zeros(2,2);
  out(i).mu = X;

  %field of view contains the point
  contains = true;
  for [val, key] = cam(i).normal_vectors
    if (X-cam(i).T)' * val <= 0
      contains = false;
      break;
    end
  end
  if ~contains
    continue;
  end


  %camera and point distance
  t = norm(X-cam(i).T);

  RT = [ cam(i).R  cam(i).T ;
           0   0      1     ] ;

  X_w_h = [ X ;
            1 ] ;

  X_c_h = RT * X_w_h;



  alfa = atan( X_c_h(2) / X_c_h(1) );

  sig = t * cos(alfa)^2 / cam(i).f_mm * cam(i).e_mm;

  Ci = [ 0     0    ;
         0  1/sig^2 ] ;

  %world <= camera <= Ci
  Rot = cam(i).R * Rot2D(alfa)';

  out(i).Ci = Rot * Ci * Rot';
  out(i).mu = X;

end

