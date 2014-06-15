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

b1 = (pi + a2 - 2*a1) / 3;
b2 = (pi + a1 - 2*a2) / 3;


figure
hold on
ezsurf(f, [0 pi/2 0 pi/2])
title('f')
plot(b1, b2, 'k*')
%axis equal
hold off


figure
hold on
ezsurf(diff(f,g1), [0 pi/2 0 pi/2])
plot(b1, b2, 'k*')
title('df / dg1')
%axis equal
hold off


figure
hold on
ezsurf(diff(f,g2), [0 pi/2 0 pi/2])
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

