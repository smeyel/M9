function [] = main10()

% M1b measurement evaluation

close all
clear
clc


myAddPath

version = 1;
xu = getX('M1b_out', version);
xd = getX('M1b_out_dist', version);
x = xu;

width_eff = max(x)-min(x);
width = 640;


% raw data
figure
plot(x);

% moving average
count = numel(x);
window= round(count / width_eff / 2);
simple = tsmovavg(x,'s',window,1);
simple = round(simple);
half=floor(window/2);
simple=[simple(half+1:end) ; simple(1:half)]; % shift with a half window
for ii = 1:numel(simple)
    simple(ii) = min(simple(ii:end));
end

figure
hold on
plot(x, 'b');
plot(simple, 'r');
hold off


% step widths
steps = zeros(width,1);
for ii = 1:numel(steps)
    steps(ii) = numel(simple(simple==ii));
end
figure
plot(steps);



function x = getX(from, version)

data = getData(from);
v = data.valid{version};
x = data.x{version};
y = data.y{version};

index_ok = 235<y & y<245 & v==1;
x = x(index_ok);

[maxx maxi] = max(x);
x = x(1:maxi);
