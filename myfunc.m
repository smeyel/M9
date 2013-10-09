
function f = myfunc(x)

E = 10;
F = 90;
Gr = 0;
Hr = 1000;

t2 = x(1)^2 + x(2)^2;
K4 = (E-F)*(Gr-Hr) + Gr*Hr;
K2 = E*Hr + F*Gr;
K0 = E*F;

nW = x(2)^2 / t2^2 * K4 + ...
     1 / t2 * K2 + ...
     K0;
f = -nW;



