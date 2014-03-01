function [] = main10()

% M1b measurement evaluation

close all
clear
clc


myAddPath

version = 1;
xu = getX('M1b_out', version);
xd = getX('M1b_out_dist', version);

width = 640;


% raw data
figure
plot(xu);
figure
plot(xd);

% moving average
count = numel(xu);
window = floor(count / width);
simple = tsmovavg(xu,'s',window,1);
half=floor(window/2);
simple=[simple(half+1:end) ; simple(1:half)]; % shift with a half window
figure
hold on
plot(simple, 'g');
plot(xu, 'b');
hold off

% raw data differences from the moving average
figure
plot(simple-xu);



function x = getX(from, version)

data = getData(from);
v = data.valid{version};
x = data.x{version};
y = data.y{version};

index_ok = 235<y & y<245 & v==1;
x = x(index_ok);

[maxx maxi] = max(x);
x = x(1:maxi);
