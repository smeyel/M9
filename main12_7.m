function [] = main12_7()

% two cameras are placed on two parametric lines
% parametric calculation

%close all
clear
clc


syms k x0 y0 positive

p21 = [x0;0];
p22 = [0;-y0];
v1 = [1;0];
v2 = (p21-p22);
on_1 = 'y==k';
on_2 = 'y==y0/x0*x-y0';
on_shift1 = 'y==0';
on_shift2 = 'y==y0/x0*x';
norm_1 = 'dot(v1,[x;y])==0';
norm_2 = 'dot(v2,[x;y])==0';


syms x y real

% t1
S = solve(eval(norm_1), eval(on_1), x, y);
t1 = [S.x ; S.y];

% t2
S = solve(eval(norm_2), eval(on_2), x, y);
t2 = [S.x ; S.y];

% m
S = solve(eval(on_1), eval(on_2), x, y);
m = [S.x ; S.y];

% m1
S = solve(eval(on_1), eval(on_shift2), x, y);
m1 = [S.x ; S.y];

% m2
S = solve(eval(on_2), eval(on_shift1), x, y);
m2 = [S.x ; S.y];

l = sqrt(sum(t2.^2));
k_ = sqrt(sum(m1.^2));
l_ = sqrt(sum(m2.^2));

% 90 deg and equal distance
% it can't be done if:
% - l_<k
% - k_<l


syms a1 a2 positive
v1 = m1-t1;
v2 = m2-t2;
p1 = t1+a1*v1;
p2 = t2+a2*v2;
eq_length = 'sqrt(sum(p1.^2)) == sqrt(sum(p2.^2))';
eq_90 = 'dot(p1,p2)==0';
S = solve(eval(eq_90), eval(eq_length), a1, a2);
p1 = subs(p1, {a1,a2}, {S.a1,S.a2});
p2 = subs(p2, {a1,a2}, {S.a1,S.a2});



syms d1 d2 real
syms a positive
v1 = m1-t1;
v2 = m2-t2;
v1 = v1 / sqrt(sum(v1.^2));
v2 = v2 / sqrt(sum(v2.^2));
pd1 = p1 + d1*v1;
pd2 = p2 + d2*v2;
pv1 = p1 + a*(pd1-p1);
pv2 = p2 + a*(pd2-p2);

Ci1 = calc_Ci(pv1);
Ci2 = calc_Ci(pv2);
Ci = Ci1 + Ci2;
return

C = inv(Ci);
lambda = eig(C);


return

par = {200,300,300};

subs(k, {k,x0,y0}, par);
subs(l, {k,x0,y0}, par);
subs(k_, {k,x0,y0}, par);
subs(l_, {k,x0,y0}, par);

subs(p1, {k,x0,y0}, par) .* [1 ; -1];
subs(p2, {k,x0,y0}, par) .* [1 ; -1];


return


function Ci = calc_Ci(p)
szog = atan(p(2)/p(1));
tav = 1/sum(p.^2);
r3 = Rot3D('z', szog);
r2 = r3(1:2,1:2);
matr = diag([0 tav]);
Ci = r2'*matr*r2;



