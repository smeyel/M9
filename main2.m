
%file for trying the actual problem

clear
clc

%warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath


[nX,nY] = meshgrid(-49:2:49, -49:2:49);


 
 
nW = (-1) * arrayfun(@(nx,ny) myfunc([nx;ny]), nX, nY);



%%3D
figure(1); clf;
contour(nX,nY,nW,900:10:1100);
axis('equal');
xlabel('x');
ylabel('y', 'rotation', 0)






minX = -40;
maxX = -15;
minY = -30;
maxY = 30;
startX = -20;
startY = -20;

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@myfunc, ...
    [startX;startY], ...
    [], [], ...
    [], [], ...
    [minX;minY], [maxX;maxY]);

oArea = [minX minY ; ...
         minX maxY ; ...
         maxX maxY ; ...
         maxX minY];

figure(1)
hold on
drawPolygon(oArea)
plot(x(1), x(2), 'r*')
hold off




%save
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), 'figures/contour_add_one_camera.eps')

