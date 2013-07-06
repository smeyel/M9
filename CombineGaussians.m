function out = CombineGaussians(in)

%in: nx1 struct array containing the fields: Ci, mu
%out: struct containing the fields: C, Ci, mu

%Ci3 = Ci1 + Ci2
%C3 = inv(Ci3)
%mu3 = C3 * (Ci1 * mu1 + Ci2 * mu2)

out.Ci = zeros();
out.mu = zeros();

for i = 1:length(in);
  out.Ci += in(i).Ci;
  out.mu += in(i).Ci * in(i).mu;
end
out.C = inv(out.Ci);
out.mu = out.C * out.mu;

