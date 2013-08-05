function out = CalculateCovariance(cam, X)

%cam: nx1 struct array of camera structs
%X: 2x1 vector
%out: nx1 struct array containing the fields: Ci, mu

global useFoV;


for i = 1:length(cam)

  out(i).Ci = zeros(2,2);
  out(i).mu = X;

  %field of view contains the point
  if useFoV
    contains = true;
    for [val, key] = cam(i).normal_vectors
      if (X-cam(i).pos)' * val <= 0
        contains = false;
        break;
      end
    end
    if ~contains
      continue;
    end
  end


  %camera and point distance
  t = norm(X-cam(i).pos);

  X_w_h = [ X ;
            1 ] ;

  X_c_h = cam(i).RT * X_w_h;
  x = X_c_h(1);
  y = X_c_h(2);

  alfa = GetAlpha2D(x,y);

  sig = t * cos(alfa)^2 / cam(i).f_mm * cam(i).e_mm;

  Ci = [ 0     0    ;
         0  1/sig^2 ] ;

  %Ci <= camera <= world
  Rot = Rot2D(-alfa) * cam(i).R;

  out(i).Ci = Rot' * Ci * Rot;
  out(i).mu = X;

end

