function alpha = GetAlpha2D(dx, dy)

% dx, dy: delta x and delta y
% alpha: 0..2*pi

alpha = 0;
if dx==0
  if dy>0
    alpha = pi/2;
  elseif dy<0
    alpha = -pi/2;
  end
else
  alpha = atan(dy/dx);
  if dx<0
    alpha = alpha + pi;
  end
end

if alpha<0
  alpha = alpha + 2*pi;
end

