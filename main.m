
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


%the columns of X are the position vectors
%X = [75 60 65 ;
%      6  0  0 ] ;

[gX,gY] = meshgrid(40:2:80, -10:2:10);
gXs = prod(size(gX));
gYs = prod(size(gY));
X = [ reshape(gX, 1, gXs) ;
      reshape(gY, 1, gYs) ] ;

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

