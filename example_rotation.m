
close all
clear
clc

myAddPath

C = diag([1,5,20]);

figure
error_ellipse(C, zeros(3,3), 'conf', 0.95);
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');


z = -pi/4;
y = pi/4;
x = 0;

z = pi/4;
y = 0;
x = 0;
Rz = [cos(z) -sin(z) 0 ; sin(z) cos(z) 0 ; 0 0 1];
Ry = [cos(y) 0 sin(y) ; 0 1 0 ; -sin(y) 0 cos(y)];
Rx = [1 0 0 ; 0 cos(x) -sin(x) ; 0 sin(x) cos(x)];
C2 = Rz*Ry*Rx*C*Rx'*Ry'*Rz';

figure
error_ellipse(C2, zeros(3,3), 'conf', 0.95);
axis('equal');
xlabel('x');
ylabel('y');
zlabel('z');