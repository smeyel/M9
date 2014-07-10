function [] = main12_7p()

% two cameras are placed on two parametric parallel lines
% parametric calculation

close all
clear
clc


eig_or_det = false;
between_or_outside = false;

m1 = 2;
m2 = 3;

syms g1 g2 real

d1 = m1/sin(g1);
d2 = m2/sin(g2);

d12 = d1^2;
d22 = d2^2;

if eig_or_det
    b1 = atan(m1/m2);
    b2 = atan(m2/m1);
    if between_or_outside
        f = 1/d12 + 1/d22 - sqrt(1/d12^2 + 1/d22^2 + 2/(d12*d22) * cos(2*(g1+g2)));
    else
        f = 1/d12 + 1/d22 - sqrt(1/d12^2 + 1/d22^2 + 2/(d12*d22) * cos(2*(g1-g2)));
        b2 = pi - b2;
    end
else
    b = linsolve([2,1;1,2],[pi;pi]);
    b1 = b(1);
    b2 = b(2);
    if between_or_outside
        f = sin(g1+g2)^2 / (d12 * d22);
    else
        f = sin(g1-g2)^2 / (d12 * d22);
        b2 = pi - b2;
    end
end


figure
hold on
ezcontour(f, [0 pi 0 pi])
plot(b1, b2, 'k*')
title('f')
%axis equal
hold off

return


