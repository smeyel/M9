
%file for trying the actual problem

%draw 2 cameras and the common covariance ellipse

clear
clc

warning('off', 'Octave:possible-matlab-short-circuit-operator');

myAddPath

figure(1); clf;
hold on;
axis([-10 200 -100 100], "equal");

%%--- camera ---
camsetup = 2;

if camsetup==1
  cam(1) = CreateCamera(0, [10;0]);
  cam(2) = CreateCamera(pi, [80;0]);
elseif camsetup==2
  cam(1) = CreateCamera(-pi/4, [10;50]);
  cam(2) = CreateCamera(pi/4, [10;-50]);
end

DrawCamera(cam)


%%--- points ---
%the columns of X are the position vectors
pointsetup = 2;

if pointsetup==1
  X = [ 75 60 65 ;
         6  0  0 ] ;
elseif pointsetup==2
  [gX,gY] = meshgrid(10:5:80, -20:2:20);
  gXs = prod(size(gX));
  gYs = prod(size(gY));
  X = [ reshape(gX, 1, gXs) ;
        reshape(gY, 1, gYs) ] ;
end



for i=1:size(X,2)
  Xt = X(:,i);
  Ci_mu = CalculateCovariance(cam, Xt);
  Comb = CombineGaussians(Ci_mu);
  Ct = Comb.C;
  Cit = Comb.Ci;
  if all(eig(Cit) > 1e-10) && all(eig(Ct) > 1e-10)
    h = my_2D_error_ellipse(Ct, Xt, 'conf', 0.95);
  end
end

hold off

