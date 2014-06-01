close all
clear
clc



a1 = 24*pi/180;
a2 = 20*pi/180;

m = 1;


syms g1 g2 real
%syms a1 a2 real

t1 = tan(a1);
t2 = tan(a2);
tgb = t1/t2 * (t2-1)/(t1-1);
b = atan(tgb);

b1 = b;
b2 = pi/2-b;

d1 = m*sin(a1)/sin(a1+g1);
d2 = m*sin(a2)/sin(a2+g2);

d1c = m*sin(a1)/sin(a1+b1);
d2c = m*sin(a2)/sin(a2+b2);

D = d1*d2/(d1c*d2c);

A = g1+g2;

K = 1/d1^4 + 1/d2^4 + 2/(d1^2*d2^2)*cos(2*A);

f = 1/d1^2 + 1/d2^2 - sqrt(K);

min1 = min(pi/2-a1, b1);
max1 = max(pi/2-a1, b1);
min2 = min(pi/2-a2, b2);
max2 = max(pi/2-a2, b2);


func = 2*D - D^2 - sin(g1+g2)^2;

func = 1/tan(a1+g1)*(D-1)-sin(2*(g1+g2));


ezsurf(func, [min1 max1 min2 max2])
title('func')
return
