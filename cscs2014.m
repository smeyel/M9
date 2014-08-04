%function [] = cscs2014()

% todo


%% preparations

close all
clear
clc

myAddPath
addpath('old/');


cams{1} = CreateCamera('pos', [0;-400]);
cams{2} = CreateCamera('pos', [0;400]);

K = cams{1}.fx / cams{1}.e;



%% merge

objX = [700;150];
newSig = 0.015;

Ci1 = calc_Ciw(cams(1), objX);
[V D] = eig(Ci1);
D(1,1) = newSig;
Ci1 = V*D*V';
C1 = pinv(Ci1);

Ci2 = calc_Ciw(cams(2), objX);
[V D] = eig(Ci2);
D(1,1) = newSig;
Ci2 = V*D*V';
C2 = pinv(Ci2);

Ci = Ci1 + Ci2;
C = pinv(Ci);

fig_merge = figure;
hold on
plot(objX(1), objX(2), 'kx')
DrawCamera(cams{1}, 'k');
DrawCamera(cams{2}, 'k');

my_2D_error_ellipse(1000*C,  objX, 'conf', 0.95, 'style', 'k');
my_2D_error_ellipse(1000*C1,  objX, 'conf', 0.95, 'style', 'k');
my_2D_error_ellipse(1000*C2,  objX, 'conf', 0.95, 'style', 'k');

axis equal
xlim([-100 1300])
ylim([-550 550])
xlabel('x')
ylabel('y', 'rotation', 0)
hold off

font = 14;
set(findall(fig_merge,'type','text'), 'FontSize', font);
set(findall(fig_merge,'type','axes'), 'FontSize', font);

export_fig(fig_merge, 'figures/cscs2014/merge.pdf', '-transparent', '-gray');



%% ellipses

[gX,gY] = meshgrid(100:100:1200, -500:100:500);

fig_ellipses = figure;
hold on

arrayfun(@(c) DrawCamera(c, 'k'), cams);
arrayfun(@(gx, gy) my_2D_error_ellipse(30*pinv(calc_Ciw(cams, [gx;gy])), [gx;gy], 'conf', 0.95, 'style', 'k'), ...
    gX, gY);

axis equal
xlim([-100 1300])
ylim([-550 550])
xlabel('x')
ylabel('y', 'rotation', 0)
hold off

font = 14;
set(findall(fig_ellipses,'type','text'), 'FontSize', font);
set(findall(fig_ellipses,'type','axes'), 'FontSize', font);

export_fig(fig_ellipses, 'figures/cscs2014/ellipses.pdf', '-transparent', '-gray');



%% wellness plots

syms x y real

x1 = cams{1}.pos(1) - x;
x2 = cams{2}.pos(1) - x;
y1 = cams{1}.pos(2) - y;
y2 = cams{2}.pos(2) - y;

d1 = sqrt(x1^2+y1^2);
d2 = sqrt(x2^2+y2^2);
g1=atan(y1/x1);
g2=atan(y2/x2);

fe = 1/d1^2 + 1/d2^2 - sqrt(1/d1^4 + 1/d2^4 + 2/(d1^2*d2^2)*cos(2*(g1-g2)));
fd = 1/d1^2 * 1/d2^2 * sin(g1-g2)^2;
ft = 1/d1^2 + 1/d2^2;

fe = K^2 * fe;
fd = K^2 * fd;
ft = K^2 * ft;

fig_fe = figure;
ezsurf(fe, [  0 1000 -1000 1000])
xlim([0 1000])
ylim([-1000 1000])
zlim([0 5])
xlabel('x')
ylabel('y', 'rotation', 0)
title('q_{eig}')

fig_fd = figure;
ezsurf(fd, [100 1000 -1000 1000])
xlim([0 1000])
ylim([-1000 1000])
zlim([0 1.5e-4])
xlabel('x')
ylabel('y', 'rotation', 0)
title('q_{det}')

fig_ft = figure;
ezsurf(ft, [100 1000 -1000 1000])
xlim([0 1000])
ylim([-1000 1000])
zlim([0 90])
xlabel('x')
ylabel('y', 'rotation', 0)
title('q_{trace}')

font = 20;
set(findall(fig_fe,'type','text'), 'FontSize', font);
set(findall(fig_fe,'type','axes'), 'FontSize', font);
set(findall(fig_fd,'type','text'), 'FontSize', font);
set(findall(fig_fd,'type','axes'), 'FontSize', font);
set(findall(fig_ft,'type','text'), 'FontSize', font);
set(findall(fig_ft,'type','axes'), 'FontSize', font);

export_fig(fig_fe, 'figures/cscs2014/fe.pdf', '-transparent', '-gray');
export_fig(fig_fd, 'figures/cscs2014/fd.pdf', '-transparent', '-gray');
export_fig(fig_ft, 'figures/cscs2014/ft.pdf', '-transparent', '-gray');



%% placement

syms g d real
R_gen = [cos(g) -sin(g) ; sin(g) cos(g)];
Ci_gen = R_gen * [1/d^2 0 ; 0 0] * R_gen';


a1 = 14*pi/180;
a2 = 10*pi/180;

m = 1;

syms g1 g2 real

d1 = m*sin(a1)/sin(a1+g1);
d2 = m*sin(a2)/sin(a2-g2);

Ci1 = subs(Ci_gen, {g,d}, {g1,d1});
Ci2 = subs(Ci_gen, {g,d}, {g2,d2});
Ci = Ci1+Ci2;


fe = 1/d1^2 + 1/d2^2 - sqrt(1/d1^4 + 1/d2^4 + 2/(d1^2*d2^2)*cos(2*(g1-g2)));
b1e =  atan((1-1/tan(a2)) / (1-1/tan(a1)));
b2e = -atan((1-1/tan(a1)) / (1-1/tan(a2)));

fd = det(Ci);
bd = (1/3 * [2 -1 ; -1 2] * [pi-a1 ; pi-a2]);
b1d =  bd(1);
b2d = -bd(2);

ft = trace(Ci);
b1t =  (pi/2 - a1);
b2t = -(pi/2 - a2);


[x y] = pol2cart(b1e, subs(d1,g1,b1e));
p1e = [x;y];
[x y] = pol2cart(b2e, subs(d2,g2,b2e));
p2e = [x;y];

[x y] = pol2cart(b1d, subs(d1,g1,b1d));
p1d = [x;y];
[x y] = pol2cart(b2d, subs(d2,g2,b2d));
p2d = [x;y];

[x y] = pol2cart(b1t, subs(d1,g1,b1t));
p1t = [x;y];
[x y] = pol2cart(b2t, subs(d2,g2,b2t));
p2t = [x;y];


fig_placement = figure;
hold on
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')

plot([0 2*m],  m*tan(a1)*[1 -1], 'k')
plot([0 2*m], -m*tan(a2)*[1 -1], 'k')

plot([0 p1e(1)], [0 p1e(2)], 'k--')
plot([0 p2e(1)], [0 p2e(2)], 'k--')
plot([0 p1d(1)], [0 p1d(2)], 'k-.')
plot([0 p2d(1)], [0 p2d(2)], 'k-.')
plot([0 p1t(1)], [0 p1t(2)], 'k:')
plot([0 p2t(1)], [0 p2t(2)], 'k:')

h1 = plot(p1e(1),p1e(2), 'ks', 'MarkerFaceColor','w');
     plot(p2e(1),p2e(2), 'ks', 'MarkerFaceColor','w')
h2 = plot(p1d(1),p1d(2), 'ko', 'MarkerFaceColor','w');
     plot(p2d(1),p2d(2), 'ko', 'MarkerFaceColor','w')
h3 = plot(p1t(1),p1t(2), 'kv', 'MarkerFaceColor','w');
     plot(p2t(1),p2t(2), 'kv', 'MarkerFaceColor','w')

axis equal
xlim([-0.1 1.1])
ylim([-0.2 0.3])
xlabel('x')
ylabel('y', 'rotation', 0)
legend([h1 h2 h3], {'eig', 'det', 'trace'})
legend boxoff
hold off

export_fig(fig_placement, 'figures/cscs2014/placement.pdf', '-transparent', '-gray');



%% objective function plots
limx = [-a1 ; pi-a1];
limy = [-(-a2) ; -(pi-a2)];
min1 = min(limx);
max1 = max(limx);
min2 = min(limy);
max2 = max(limy);


fig_cfe = figure;
hold on
ezcontour(fe, [min1 max1 min2 max2])
h = ezplot('g1-g2=0'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=pi/2'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=pi'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=3*pi/2'); set(h,'LineStyle', ':');
title('q_{eig}')
plot(b1e, b2e, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
axis equal
xlim([min1 max1])
ylim([min2 max2])
xlabel('$\alpha_1$','interpreter','latex')
ylabel('$\alpha_2$','interpreter','latex', 'rotation', 0)
hold off


fig_cfd = figure;
hold on
ezcontour(fd, [min1 max1 min2 max2])
h = ezplot('g1-g2=0'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=pi/2'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=pi'); set(h,'LineStyle', ':');
h = ezplot('g1-g2=3*pi/2'); set(h,'LineStyle', ':');
title('q_{det}')
plot(b1d, b2d, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
axis equal
xlim([min1 max1])
ylim([min2 max2])
xlabel('$\alpha_1$','interpreter','latex')
ylabel('$\alpha_2$','interpreter','latex', 'rotation', 0)
hold off


fig_cft = figure;
hold on
ezcontour(ft, [min1 max1 min2 max2])
title('q_{trace}')
plot(b1t, b2t, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
axis equal
xlim([min1 max1])
ylim([min2 max2])
xlabel('$\alpha_1$','interpreter','latex')
ylabel('$\alpha_2$','interpreter','latex', 'rotation', 0)
hold off

font = 18;
set(findall(fig_cfe,'type','text'), 'FontSize', font);
set(findall(fig_cfe,'type','axes'), 'FontSize', font);
set(findall(fig_cfd,'type','text'), 'FontSize', font);
set(findall(fig_cfd,'type','axes'), 'FontSize', font);
set(findall(fig_cft,'type','text'), 'FontSize', font);
set(findall(fig_cft,'type','axes'), 'FontSize', font);

export_fig(fig_cfe, 'figures/cscs2014/cfe.pdf', '-transparent', '-gray');
export_fig(fig_cfd, 'figures/cscs2014/cfd.pdf', '-transparent', '-gray');
export_fig(fig_cft, 'figures/cscs2014/cft.pdf', '-transparent', '-gray');

