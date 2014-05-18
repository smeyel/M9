function [] = main12_6()

% two cameras are placed on two parametric lines
% surf is plotted

close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath


global segment
segment = {};
%segment{1} = [-500 -500 -500  500];
%segment{2} = [ 500 -500  500  500];
%segment{1} = [-1000   500 1000  500];
%segment{2} = [ -500 -1000 -500 1000];
%segment{1} = [ -400  -100    0 -100];
%segment{2} = [ -400     0 -300  100];
k = -200;
segment{1} = [0 k -k k];
segment{2} = [150 150 300 0];


x0 = cellfun(@(s) s(1), segment)';
y0 = cellfun(@(s) s(2), segment)';
x1 = cellfun(@(s) s(3), segment)';
y1 = cellfun(@(s) s(4), segment)';
dx = x1-x0;
dy = y1-y0;

g1 = linspace(0,1,16);
g2 = linspace(0,1,16);
[t1,t2] = ndgrid(g1,g2);

tW = arrayfun(@(t1,t2) real(myfun([x0+[t1;t2].*dx,y0+[t1;t2].*dy])), t1, t2);

figure
hold on
surf(t1,t2,tW)
xlabel('t1');
ylabel('t2');
grid on
hold off


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));


