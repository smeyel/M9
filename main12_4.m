function [] = main12_4()

% measure calculated on an ndgrid defined on the main12_3 segments


close all
clear
clc

global useDetectAngle
useDetectAngle = false;

myAddPath

global segment
segment{1} = [-500 -400 400 -600];
segment{2} = [600 -300 -200 400];
segment{3} = [-700 0 200 600];
segment{4} = [-400 -400 -400 400];





x0 = cellfun(@(s) s(1), segment)';
y0 = cellfun(@(s) s(2), segment)';
x1 = cellfun(@(s) s(3), segment)';
y1 = cellfun(@(s) s(4), segment)';
dx = x1-x0;
dy = y1-y0;

count = numel(segment);

g1 = linspace(0,1,15);
g2 = linspace(0,1,15);
g3 = linspace(0,1,15);
g4 = linspace(0,1,15);
[t1,t2,t3,t4] = ndgrid(g1,g2,g3,g4);

tW = arrayfun(@(t1,t2,t3,t4) myfun([x0+[t1;t2;t3;t4].*dx,y0+[t1;t2;t3;t4].*dy]), t1, t2, t3, t4);
mtW = max(max(max(max(tW))));
save('main12_4')

mtW

return


function W = myfun(Xs)
x_cell = num2cell(Xs,2);
cams = cellfun(@(x) CreateCamera('pos', x'), x_cell, 'uni', false);
W = 10000*min(eig(calc_Ciw(cams, zeros(cams{1}.dim,1))));



