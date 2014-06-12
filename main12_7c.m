close all
clear
clc


a1 = 24*pi/180;
a2 = 20*pi/180;

t1 = tan(a1);
t2 = tan(a2);
tgb1 = t1/t2 * (t2-1)/(t1-1);
b1 = atan(tgb1);
tgb2 = 1/tgb1;
b2 = pi/2-b1;

syms g1 g2 real


A1 = tan(a1);
A2 = tan(a2);
K1 = 1/A1^2;
K2 = 1/A2^2;
d1r2 = K1*(A1+g1)^2/(1+g1^2);
d2r2 = K2*(A2+g2)^2/(1+g2^2);
c12 = ((g1*g2-1)^2-(g1+g2)^2) / ((1+g1^2)*(1+g2^2));
f = d1r2 + d2r2 - sqrt(d1r2^2 + d2r2^2 + 2*d1r2*d2r2*c12);

mm = 5;

figure
hold on
ezsurf(f, [0  mm 0 mm])
plot(tgb1, tgb2, 'k*')
hold off

return

