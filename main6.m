function [] = main6()

% the standard deviation of marker location estimations,
% shown for the 3 axes separately
% the cameras and measurement setup is defined
% in the ICSSE_2013 publication

close all
clear
clc

global useDetectAngle
useDetectAngle = true;


myAddPath

data = getData('M1R1_c2w');
dataNames = fieldnames(data);
for i = 1:numel(dataNames)
    cams{i} = CreateCamera('c2w', data.(dataNames{i}));
end

data = getData('M1R1_stats');
means = data.LocationMean3Ray;

cc = max(size(means)); % column count
TheoryLocationStd = zeros(cc,3);
TheoryLocationMeasure = zeros(cc);
for i=1:cc

    % observed point
    p = means(i,:)';

    if any(isnan(p))
        sigres = [0;0;0];
    else
        % covariance matrix inverses and the resulting covariace matrix
        Ciw = calc_Ciw(cams, p);
        Cw = inv(Ciw);

        % variance in the x-y-z directions
        [V D] = eig(Cw);
        sigres = V.^2 * diag(D);
        measure = max(diag(D));
    end
    TheoryLocationStd(i,1:3) = sigres';
    TheoryLocationMeasure(i) = measure;

end

figure
bar(TheoryLocationStd(:,:));
axis([1 size(TheoryLocationStd,1) 0 5]);
legend('X','Y','Z');
xlabel('Location index');
ylabel('Standard deviation');
daspect([3 1 1])
set(gca, 'ytick', 0:0.5:5)

figure
bar(TheoryLocationMeasure(:));
axis([1 size(TheoryLocationMeasure,1) 0 1.5]);
legend('measure');
xlabel('Location index');
ylabel('Measure');

