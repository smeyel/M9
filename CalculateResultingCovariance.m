function out = CalculateResultingCovariance(cam, X)

%cam: nx1 struct array of camera structs
%X: 2x1 vector
%out: struct containing the fields: valid, C, Ci, mu

out = CombineGaussians(CalculateCovariance(cam, X));

