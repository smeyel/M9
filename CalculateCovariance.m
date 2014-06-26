function out = CalculateCovariance(cam, X)

%cam: nx1 struct array of camera structs
%X: 2x1 vector
%out: nx1 struct array containing the fields: Ci, mu

global useFoV;


for i = 1:length(cam)

  out(i).Ci = zeros(2,2);
  out(i).mu = X;


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

  %cinti2013 conference article: two cameras, one observed point, figure
  Ci_kamu = [ 0.01/sig^2     0    ;
                    0     1/sig^2 ] ;

  %Ci <= camera <= world
  Rot = Rot2D(-alfa) * cam(i).R;

  out(i).Ci = Rot' * Ci * Rot;
  out(i).Ci_kamu = Rot' * Ci_kamu * Rot;
  out(i).mu = X;

end

