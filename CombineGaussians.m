function out = CombineGaussians(in)

%in: nx1 struct array containing the fields: Ci, mu
%out: struct containing the fields: valid, C, Ci, mu

%Ci3 = Ci1 + Ci2
%C3 = inv(Ci3)
%mu3 = C3 * (Ci1 * mu1 + Ci2 * mu2)

out.C = zeros();
out.Ci = zeros();
out.mu = zeros();
mmu = zeros();

for i = 1:length(in);
  out.Ci += in(i).Ci;
  mmu += in(i).Ci * in(i).mu;
end

out.valid = all(eig(out.Ci) > 1e-10);

if out.valid
  out.C = inv(out.Ci);
  out.mu = out.C * mmu;
end

