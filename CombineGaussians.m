function out = CombineGaussians(in)

%in: 1xn struct array containing the fields: Ci, mu
%out: struct containing the fields: C, Ci, mu


n = length(in);

%Ci
out.Ci = zeros();
for i=1:n
  out.Ci += in(i).Ci;
end

%C
out.C = inv(out.Ci);

%mu
out.mu = zeros();
for i=1:n
  out.mu += in(i).Ci * in(i).mu;
end
out.mu = out.C * out.mu;

