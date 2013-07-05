function DrawCamera(cam)

%cam: camera struct


default_center_circle_radius = 2;
p = circleToPolygon([cam.T(1) cam.T(2) default_center_circle_radius]);
fill(p(:,1), p(:,2), "b");

c = cam.T;
f = cam.normal_vectors.front;
dl = Rot2D(pi/2) * cam.normal_vectors.low;
if(dl' * f < 0)
  dl *= -1;
end
dh = Rot2D(pi/2) * cam.normal_vectors.high;
if(dh' * f < 0)
  dh *= -1;
end

drawRay([c;dl]');
drawRay([c;dh]');

