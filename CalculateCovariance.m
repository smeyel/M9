function out = CalculateCovariance(cam, X)

%cam: camera struct
%X: 2x1 vector
%out: struct containing the fields: Ci, mu

global useFoV;

%preallocation
out = struct('Ci', zeros(2,2), ...
             'mu', zeros(2,1));

out.Ci = zeros(2,2);
out.mu = X;

%field of view contains the point
if useFoV
  S = cam.normal_vectors;
  SNames = fieldnames(S);
  for loopIndex = 1:numel(SNames)
    val = S.(SNames{loopIndex});
    if (X-cam.pos)' * val <= 0
      return;
    end
  end
end


%camera and point distance
t = norm(X-cam.pos);

X_w_h = [ X ;
          1 ] ;

X_c_h = cam.RT * X_w_h;
x = X_c_h(1);
y = X_c_h(2);

alfa = GetAlpha2D(x,y);

sig = t * cos(alfa)^2 / cam.f_mm * cam.e_mm;

Ci = [ 0     0    ;
       0  1/sig^2 ] ;

%cinti2013 conference article: two cameras, one observed point, figure
Ci_kamu = [ 0.01/sig^2     0    ;
                  0     1/sig^2 ] ;

%Ci <= camera <= world
Rot = Rot2D(-alfa) * cam.R;

out.Ci = Rot' * Ci * Rot;
out(i).Ci_kamu = Rot' * Ci_kamu * Rot;
out.mu = X;

