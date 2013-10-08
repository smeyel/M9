function out = CalculateResultingCovariance(cams, X)

%cams: nx1 struct array of camera structs
%X: 2x1 vector
%out: struct containing the fields: valid, C, Ci, mu

covs = arrayfun(@(cam) CalculateCovariance(cam, X), cams);
out = CombineGaussians(covs);

