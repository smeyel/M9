function data = getData(from)

if strcmp(from, 'M1R1_c2w')
    % from the ICSSE_2013 publication
    % T matrices from the cpp program print
    % transformation: camera => world
    data.c2w0 = [-0.69097477,  -0.13459364,  0.71023828, -674.35431;
                  0.011185364, -0.98438656, -0.17566407,  263.93491;
                  0.72279227,  -0.11343517,  0.68169177, -677.69617;
                  0,            0,           0,             1];
    data.c2w1 = [-0.70978242,  -0.21733436, -0.67005581,  871.03137;
                 -0.080154344, -0.92011875,  0.38334945, -266.29767;
                 -0.69984591,   0.32580256,  0.63566369, -522.82904;
                  0,            0,           0,             1];
    data.c2w2 = [-0.97232783,   0.11579948, 0.20290148,  -28.091051;
                 -0.020656202, -0.90772098, 0.41906551, -265.4436;
                  0.23270552,   0.40327793, 0.88499433, -758.5166;
                  0,            0,          0,             1];
elseif strcmp(from, 'M1R1_stats')
    indata = dlmread('M1R1_stableframeDataProcessor_stats.csv',' ');
    % Statistics for every stable frame interval (location)
    data.LocationMean2Ray = indata(:,1:3);
    data.LocationMean3Ray = indata(:,4:6);
    data.LocationMeanAll = indata(:,7:9);
    data.LocationStd2Ray = indata(:,10:12);
    data.LocationStd3Ray = indata(:,13:15);
    data.LocationStdAll = indata(:,16:18);
    data.LocationEffectiveStd2Ray = indata(:,19);
    data.LocationEffectiveStd3Ray = indata(:,20);
    data.LocationEffectiveStdAll = indata(:,21);
elseif strcmp(from, 'ps3eye_intrinsics_red.xml')
    % image_Width
    data.width = 640;
    % image_Height
    data.height = 480;
    % Avg_Reprojection_Error
    data.e = 0.7758;
    % Camera_Matrix
    data.fx = 789.1510;
    data.fy = 789.1510;
    data.cx = 319.5;
    data.cy = 239.5;
    % Distortion_Coefficients
    data.k1 = -0.1319227908065607;
    data.k2 = 0.1525229491386591;
    data.p1 = 0;
    data.p2 = 0;
    data.k3 = -1.115072806560045;
else
    error('No such from!!!')
end



