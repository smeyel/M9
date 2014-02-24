function [] = main6()

% the standard deviation of marker location estimations,
% shown for the 3 axes separately
% the cameras and measurement setup is defined
% in the ICSSE_2013 publication

close all
clear
clc

global useDetectAngle
useDetectAngle = false;


myAddPath

data = getData('M1R1_c2w');
dataNames = fieldnames(data);
for i = 1:numel(dataNames)
    cams{i} = CreateCamera('c2w', data.(dataNames{i}));
end

data = getData('M1R1_stats');
means = data.LocationMean3Ray;
stdall = data.LocationStdAll;

cc = max(size(means)); % column count
TheoryLocationStd = zeros(cc,3);
TheoryLocationMeasure = zeros(cc);
for i=1:cc

    % observed point
    p = means(i,:)';

    if any(isnan(p))
        sigres = [NaN;NaN;NaN];
    else
        % covariance matrix inverses and the resulting covariace matrix
        Ciw = calc_Ciw(cams, p);
        Cw = inv(Ciw);

        % deviation in the w direction:
        % w: column vector, the requested direction
        % v(i): column vectors, the eigenvectors of the covariance matrix
        % d(i): number, the variances in the v(i) directions
        % V: matrix, the v(i) vectors are its column vectors
        % D: matrix, diagonal matrix with d(i) elements
        % sig2w: number, the variation in the w direction
        % sigw: number, the deviation in the w direction
        % sig2w = sum( d(i)^2 * (v(i)*w)^2 ) =
        %       = [d(1);d(2);d(3)]'.^2 * ([v(1);v(2);v(3)]'*w).^2 =
        %       = (w'*[v(1);v(2);v(3)]).^2 * [d(1);d(2);d(3)].^2 =
        %       = (w'*V).^2 * diag(D)
        % sigw = sqrt(sig2w) = sqrt( (w'*V).^2 * diag(D) )
        % diag(D) is a column vector and it multiplies one row vector
        % the result is a number
        %
        % deviation in the w1-w2-w3 directions:
        % W: matrix, the wi vectors are its column vectors
        % sig2W: column vector, the variations in the w1-w2-w3 direction
        % sigW: column vector, the deviations in the w1-w2-w3 direction
        % sig2W = ([w1;w2;w3]'*V).^2 * diag(D) =
        %       = (W'*V).^2 * diag(D)
        % sigW = sqrt(sig2W) = sqrt( (W'*V).^2 * diag(D) )
        %
        % deviation in the x-y-z directions
        % W = eye(3)
        [V D] = eig(Cw);
        sigres = sqrt( V.^2 * diag(D) );
        measure = max(diag(D));
    end
    TheoryLocationStd(i,1:3) = sigres';
    TheoryLocationMeasure(i) = measure;

end

figure
bar(stdall(:,[1 3 2]));
axis([1 size(stdall,1) 0 5]);
legend('X','Y','Z');
xlabel('Location index');
ylabel('Standard deviation');
title('Standard deviation (measurement)');
daspect([3 1 1])
set(gca, 'ytick', 0:0.5:5)

figure
bar(TheoryLocationStd(:,:));
axis([1 size(TheoryLocationStd,1) 0 5]);
legend('X','Y','Z');
xlabel('Location index');
ylabel('Standard deviation');
title('Standard deviation (theory)');
daspect([3 1 1])
set(gca, 'ytick', 0:0.5:5)

diffstd = stdall(:,[1 3 2]) - TheoryLocationStd;
figure
bar(diffstd(:,:));
axis([1 size(diffstd,1) -1 1.5]);
legend('X','Y','Z');
xlabel('Location index');
ylabel('Standard deviation');
title('Standard deviation (measurement-theory)');
daspect([6 1 1])
set(gca, 'ytick', -1:0.5:1.5)

figure
bar(TheoryLocationMeasure(:));
axis([1 size(TheoryLocationMeasure,1) 0 1.5]);
legend('measure');
xlabel('Location index');
ylabel('Measure');
title('Measure');

