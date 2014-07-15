function [] = main13()

% two cameras are placed on two parametric lines
% parametric calculation
% determinant measure

close all
clear
clc



syms g d real
R_gen = [cos(g) -sin(g) ; sin(g) cos(g)];
Ci_gen = R_gen * [0 0 ; 0 1/d^2] * R_gen';



% general case
if false
syms g1 d1 real
syms g2 d2 real
syms g3 d3 real

Ci1 = subs(Ci_gen, {g,d}, {g1,d1});
Ci2 = subs(Ci_gen, {g,d}, {g2,d2});
Ci3 = subs(Ci_gen, {g,d}, {g3,d3});

Ci_res = Ci1+Ci2+Ci3;
f = det(Ci_res);
pretty(f)
return
end



a1 = 24*pi/180;
a2 = 20*pi/180;

m = 1;


syms g1 g2 real

d1 = m*sin(a1)/sin(a1+g1);
d2 = m*sin(a2)/sin(a2+g2);

Ci1 = subs(Ci_gen, {g,d}, {g1,d1});
Ci2 = subs(Ci_gen, {g,d}, {-g2,d2});
Ci = Ci1+Ci2;

f = det(Ci);

b = (1/3 * [2 -1 ; -1 2] * [pi-a1 ; pi-a2]);
b1 = b(1);
b2 = b(2);
b1
b2

b1_e = atan((1-1/tan(a2)) / (1-1/tan(a1)));
b2_e = atan((1-1/tan(a1)) / (1-1/tan(a2)));


[x y] = pol2cart(b1_e, subs(d1,g1,b1_e));
p1_e = [x;y];

[x y] = pol2cart(-b2_e, subs(d2,g2,b2_e));
p2_e = [x;y];

[x y] = pol2cart(b1, subs(d1,g1,b1));
p1 = [x;y];

[x y] = pol2cart(-b2, subs(d2,g2,b2));
p2 = [x;y];

[x y] = pol2cart(pi/2-a1, subs(d1,g1,pi/2-a1));
talp1 = [x;y];

[x y] = pol2cart(-(pi/2-a2), subs(d2,g2,pi/2-a2));
talp2 = [x;y];

[x y] = pol2cart(a2, subs(d1,g1,a2));
metszet1 = [x;y];

[x y] = pol2cart(-a1, subs(d2,g2,a1));
metszet2 = [x;y];


figure
plot([0 m],[0 0])
hold on
plot([0 m],[m*tan(a1) 0])
plot([0 m],[-m*tan(a2) 0])
plot([0 p1(1)],[0 p1(2)])
plot([0 p2(1)],[0 p2(2)])

plot([p1(1) p2(1)],[p1(2) p2(2)], 'r')
plot(p1(1),p1(2), 'r*')
plot(p2(1),p2(2), 'r*')
plot(p1_e(1),p1_e(2), 'c*')
plot(p2_e(1),p2_e(2), 'c*')
plot([0 talp1(1)], [0 talp1(2)], 'k')
plot([0 talp2(1)], [0 talp2(2)], 'k')
plot([0 metszet1(1)], [0 metszet1(2)], 'g')
plot([0 metszet2(1)], [0 metszet2(2)], 'g')

axis equal
hold off


%return

plotLimits = 2;

if plotLimits == 1
    min1 = min(pi/2-a1, a2);
    max1 = max(pi/2-a1, a2);
    min2 = min(pi/2-a2, a1);
    max2 = max(pi/2-a2, a1);
elseif plotLimits == 2
	min1 = -a1;
	max1 = pi-a1;
	min2 = -a2;
	max2 = pi-a2;
elseif plotLimits == 3
	min1 = 0;
	max1 = pi/2;
	min2 = 0;
	max2 = pi/2;
end


figure
hold on
ezsurf(f, [min1 max1 min2 max2])
title('f')
plot(b1, b2, 'k*')
%axis equal
hold off


figure
hold on
ezsurf(diff(f,g1), [min1 max1 min2 max2])
plot(b1, b2, 'k*')
title('df / dg1')
%axis equal
hold off


figure
hold on
ezsurf(diff(f,g2), [min1 max1 min2 max2])
plot(b1, b2, 'k*')
title('df / dg2')
%axis equal
hold off


figure
hold on
h1 = ezplot(2*g1+g2+a1==pi, [0 pi/2 0 pi/2]);
set(h1, 'Color', 'b')
h2 = ezplot(2*g2+g1+a2==pi, [0 pi/2 0 pi/2]);
set(h2, 'Color', 'r')
plot(b1, b2, 'k*')
legend('df / dg1 = 0', 'df / dg2 = 0')
title('df / dg = 0')
axis equal
hold off

