
%file for trying the actual problem

clear
clc


[gX,gY] = meshgrid(-49:2:49, -49:2:49);



E = 10;
F = 90;
Gr = 0;
Hr = 1000;


%preallocation
gW = zeros(size(gX,1), ...
           size(gX,2));

%for every point in grid
for i=1:size(gX,1)
  for j=1:size(gX,2)

    x = gX(i,j);
    y = gY(i,j);

    t2 = x^2+y^2;
    K4 = (E-F)*(Gr-Hr) + Gr*Hr;
    K2 = E*Hr + F*Gr;
    K0 = E*F;

    gW(i,j) = y^2/t2^2*K4 + 1/t2*K2 + K0;

  end
end



%%3D
figure(1); clf;
contour(gX,gY,gW,900:10:1100);
axis('equal');
xlabel('x');
ylabel('y', 'rotation', 0)



%save
colormap([zeros(63,3) ; ones(1,3)]);
saveas(figure(1), 'figures/contour_add_one_camera.eps')

