function [] = main12_7()

% two cameras are placed on two parametric lines
% parametric calculation

%close all
clear
clc


syms a1 a2 m positive

d1 = sin(a1) * m;
d2 = sin(a2) * m;
on_1 = 'y==tan(-a1)*(x-m)';
on_2 = 'y==tan(a2)*(x-m)';
on_shift1 = 'y==-tan(a1)*x';
on_shift2 = 'y==tan(a2)*x';

syms x y real

% m1
S = solve(eval(on_1), eval(on_shift2), x, y);
m1 = [S.x ; S.y];

% m2
S = solve(eval(on_2), eval(on_shift1), x, y);
m2 = [S.x ; S.y];

% t1
t1 = d1 * [sin(a1); cos(a1)];

% t2
t2 = d2 * [sin(a2); cos(a2)];


m1_ = sqrt(sum(m1.^2));
m2_ = sqrt(sum(m2.^2));
t1_ = sqrt(sum(t1.^2));
t2_ = sqrt(sum(t2.^2));

% 90 deg and equal distance
% it can't be done if:
% - m1_ < t2_
% - m2_ < t1_

syms j1 j2 positive
v1 = m1-t1;
v2 = m2-t2;
p1 = t1+j1*v1;
p2 = t2+j2*v2;
eq_length = 'sqrt(sum(p1.^2)) == sqrt(sum(p2.^2))';
eq_90 = 'dot(p1,p2)==0';



v1 = m1-t1;
v2 = m2-t2;
v1 = v1 / sqrt(sum(v1.^2));
v2 = v2 / sqrt(sum(v2.^2));
pj1 = p1 + j1*v1;
pj2 = p2 + j2*v2;

Ci1 = calc_Ci(pj1);
Ci2 = calc_Ci(pj2);
Ci = Ci1 + Ci2;


return

S = solve(eval(eq_90), eval(eq_length), j1, j2);
p1 = subs(p1, {j1,j2}, {S.j1,S.j2});
p2 = subs(p2, {j1,j2}, {S.j1,S.j2});

p1
p2



function Ci = calc_Ci(p)
szog = atan(p(2)/p(1));
tav = 1/sum(p.^2);
r3 = Rot3D('z', szog);
r2 = r3(1:2,1:2);
matr = diag([0 tav]);
Ci = r2'*matr*r2;



