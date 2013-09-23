function DrawCamera(cam, color)

%cam: camera struct

global useFoV;

default_center_circle_radius = 2;

if(nargin < 2)
    color='b';
end

p = circleToPolygon([cam.pos ; default_center_circle_radius]);
fill(p(:,1), p(:,2), color);

%center
c = cam.pos;

%normal vector for front
nf = cam.normal_vectors.front;

if useFoV

  %direction low
  dl = Rot2D(pi/2) * cam.normal_vectors.low;
  if(dl' * nf < 0)
    dl = -dl;
  end

  %direction high
  dh = Rot2D(pi/2) * cam.normal_vectors.high;
  if(dh' * nf < 0)
    dh = -dh;
  end

  drawRay([c;dl]', color);
  drawRay([c;dh]', color);

else

  l = 10; %length
  plot([c(1) c(1)+l*nf(1)], [c(2) c(2)+l*nf(2)], color)

end

