function DrawCamera(cam)

%cam: nx1 struct array of camera structs

default_center_circle_radius = 2;

for i = 1:length(cam)

  p = circleToPolygon([cam(i).T(1) cam(i).T(2) default_center_circle_radius]);
  fill(p(:,1), p(:,2), "b");

  %center
  c = cam(i).T;

  %normal vector for front
  nf = cam(i).normal_vectors.front;

  %direction low
  dl = Rot2D(pi/2) * cam(i).normal_vectors.low;
  if(dl' * nf < 0)
    dl *= -1;
  end

  %direction high
  dh = Rot2D(pi/2) * cam(i).normal_vectors.high;
  if(dh' * nf < 0)
    dh *= -1;
  end

  drawRay([c;dl]');
  drawRay([c;dh]');

end

