function main9

% the variation taking the effect of the distortion into account
% using virtual cameras with parameters defined in ICSSE_2013 publication
% the measure is the smallest eigenvalue of the inverse cov. matr.
% it has to be maximized
% it is displayed in the last figure
% previous figures contain deviation mesh plots
% the deviation has to be minimized

close all
clear
clc


cam = CreateCamera('oripos', [0;0;1], [0;0;0]);
[mapx mapy] = create_map(cam);


step = 8;
figure
hold on
for u=1:step:cam.width
    plot(mapx(:,u),1:cam.height)
end
for v=1:step:cam.height
    plot(1:cam.width, mapy(v,:))
end
axis([1 640 1 480], 'equal')
hold off


[sy sx] = size(mapx);
xm1 = zeros(sy, sx);
xp1 = zeros(sy, sx);
ym1 = zeros(sy, sx);
yp1 = zeros(sy, sx);
for i=1:sy
    for j=1:sx
        if j > 1 ; xm1(i,j) = norm([mapx(i,j)-mapx(i,j-1) mapy(i,j)-mapy(i,j-1)]) ; end
        if j < sx ; xp1(i,j) = norm([mapx(i,j)-mapx(i,j+1) mapy(i,j)-mapy(i,j+1)]) ; end
        if i > 1 ; ym1(i,j) = norm([mapx(i,j)-mapx(i-1,j) mapy(i,j)-mapy(i-1,j)]) ; end
        if i < sy ; yp1(i,j) = norm([mapx(i,j)-mapx(i+1,j) mapy(i,j)-mapy(i+1,j)]) ; end
    end
end

wx = ones(sy, sx)/2;
wx(:,1) = 1;
wx(:,end) = 1;
wy = ones(sy, sx)/2;
wy(1,:) = 1;
wy(end,:) = 1;
mult_sigx = 1 ./ (wx.*(xm1+xp1));
mult_sigy = 1 ./ (wy.*(ym1+yp1));

figure
mesh(mult_sigx)
figure
mesh(mult_sigy)

x = 1:cam.width;
px = cos(atan((x-cam.cx)/cam.fx)).^2;
px = repmat(px, cam.height, 1);
sig_newx = mult_sigx.*px;
figure
mesh(sig_newx)

y = (1:cam.height)';
py = cos(atan((y-cam.cy)/cam.fy)).^2;
py = repmat(py, 1, cam.width);
sig_newy = mult_sigy.*py;
figure
mesh(sig_newy)

meas_new = max(sig_newx, sig_newy).^(-2);
figure
mesh(meas_new)
figure
contour(meas_new)



function [mapx mapy] = create_map(cam)
ir = cam.ir;
mapx = zeros(cam.height, cam.width);
mapy = zeros(cam.height, cam.width);
for i=1:cam.height
    vx = i*ir(2) + ir(3);
    vy = i*ir(5) + ir(6);
    vw = i*ir(8) + ir(9);
    for j=1:cam.width
        x = vx/vw;
        y = vy/vw;
        [mapx(i,j) mapy(i,j)] = calc_map_for_one_pixel(cam, x, y);
        vx = vx +ir(1);
        vy = vy +ir(4);
        vw = vw +ir(7);
    end
end


function [u v] = calc_map_for_one_pixel(cam, x, y)

k1 = cam.k1;
k2 = cam.k2;
p1 = cam.p1;
p2 = cam.p2;
k3 = cam.k3;
fx = cam.fx;
fy = cam.fy;
u0 = cam.cx;
v0 = cam.cy;

x2 = x*x;
y2 = y*y;
r2 = x2 + y2;
m2xy = 2*x*y;
kr = 1 + ((k3*r2 + k2)*r2 + k1)*r2;
u = fx*(x*kr + p1*m2xy + p2*(r2 + 2*x2)) + u0;
v = fy*(y*kr + p1*(r2 + 2*y2) + p2*m2xy) + v0;
