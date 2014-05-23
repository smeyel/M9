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
b*180/pi


syms g1 g2 real

d1 = m*sin(a1)/sin(a1+g1);
d2 = m*sin(a2)/sin(a2+g2);

d12 = d1^2;
d22 = d2^2;

f = 1/d12 + 1/d22 - sqrt(1/d12^2 + 1/d22^2 + 2/(d12*d22) * cos(2*(g1+g2)));

figure
hold on
ezcontour(f, [0 pi/2 0 pi/2])
plot(b, pi/2-b, 'k*')
title('f')
axis equal
hold off

figure
hold on
ezsurf(diff(f,g1), [0 pi/2 0 pi/2])
plot(b, pi/2-b, 'k*')
title('df / dg1')
axis equal
hold off

figure
hold on
ezsurf(diff(f,g2), [0 pi/2 0 pi/2])
plot(b, pi/2-b, 'k*')
title('df / dg2')
axis equal
hold off



syms x y real
df1 =  char(subs(diff(f,g1), {g1,g2}, {x,y}));
df2 =  char(subs(diff(f,g2), {g1,g2}, {x,y}));
clear x y

y = linspace(0,pi/2,200);
xmin = 0;
xmax = pi/2;
eps = 0.001;

x = arrayfun(@(y) calc_intersection(df1, y, xmin, xmax, eps), y);
figure
hold on
plot(x,y)
plot(b, pi/2-b, 'k*')
plot([pi/2-a1 pi/2-a1], [0 pi/2], 'k')
plot([a2 a2], [0 pi/2], 'g')
plot([0 pi/2], [pi/2-a2 pi/2-a2], 'k')
plot([0 pi/2], [a1 a1], 'g')
title('df / dg1 = 0')
axis equal
hold off

y = linspace(0,pi/2,200);
x = arrayfun(@(y) calc_intersection(df2, y, xmin, xmax, eps), y);
figure
hold on
plot(x,y)
plot(b, pi/2-b, 'k*')
plot([pi/2-a1 pi/2-a1], [0 pi/2], 'k')
plot([a2 a2], [0 pi/2], 'g')
plot([0 pi/2], [pi/2-a2 pi/2-a2], 'k')
plot([0 pi/2], [a1 a1], 'g')
title('df / dg2 = 0')
axis equal
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


function intersection = calc_intersection(func_str, y, xmin, xmax, eps)

low = xmin;
high = xmax;

while (high-low) > eps
    mid = (low+high)/2;
    x = mid;
    y = y;
    val = eval(func_str);
    if val > 0
        low = mid;
    elseif val < 0
        high = mid;
    else
        break;
    end
end
intersection = mid;

