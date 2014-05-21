function [] = main12_7()

% two cameras are placed on two parametric lines
% parametric calculation
% calculate beta from angles with tangents

close all
clear
clc




a1 = 60*pi/180;
a2 = 80*pi/180;

m = 10;


syms as bs real
metszet = m * [tan(as)/(tan(as)+tan(bs)) ; tan(as)*tan(bs)/(tan(as)+tan(bs))];

talp1 = subs(metszet, {as,bs}, {a1,pi/2-a1});
talp2 = subs(metszet, {as,bs}, {a2,pi/2-a2}) .* [1;-1];
metszet1 = subs(metszet, {as,bs}, {a1,a2});
metszet2 = subs(metszet, {as,bs}, {a2,a1}) .* [1;-1];

talp1_ = sqrt(sum(talp1.^2));
talp2_ = sqrt(sum(talp2.^2));
metszet1_ = sqrt(sum(metszet1.^2));
metszet2_ = sqrt(sum(metszet2.^2));



% 90 deg and equal distance
% it can't be done if:
% - metszet1_ <= talp2_
% - metszet2_ <= talp1_

if metszet1_ <= talp2_
    b = a2;
elseif metszet2_ <= talp1_
    b = pi/2 - a1;
else
    t1 = tan(a1);
    t2 = tan(a2);
    if a1 == pi/4
        tgb = t1 * (1-t2)/(t1-t2);
    elseif a2 == pi/4
        tgb = (t2-t1) / (t2 * (t1-1));
    else
        tgb = t1/t2 * (t2-1)/(t1-1);
    end
    b = atan(tgb);
end


p1 =  subs(metszet, {as,bs}, {a1,b});
p2 =  subs(metszet, {as,bs}, {a2,pi/2-b}) .* [1;-1];



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
plot([0 talp1(1)], [0 talp1(2)], 'k')
plot([0 talp2(1)], [0 talp2(2)], 'k')
plot([0 metszet1(1)], [0 metszet1(2)], 'g')
plot([0 metszet2(1)], [0 metszet2(2)], 'g')


axis equal
axis([0 m p2(2) p1(2)])
hold off

dot(p1,p2)
norm(p1)-norm(p2)
b
pi/2-b


syms g1 g2 real

d12 = m^2 * t1^2*(1+tan(g1)^2) / (t1+tan(g1))^2;
d22 = m^2 * t2^2*(1+tan(g2)^2) / (t2+tan(g2))^2;

f = 1/d12 + 1/d22 - sqrt(1/d12^2 + 1/d22^2 + 2/(d12*d22) * cos(2*(g1+g2)));

figure
hold on
ezcontour(f, [0 pi/2 0 pi/2])
plot(b, pi/2-b, 'k*')
hold off


return




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



