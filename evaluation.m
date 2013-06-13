% Reset MATLAB environment.
clear all;
close all;

% Add algorithms.
addpath('algorithms/pixel/');

% Load data.
data = load('data/data.mat');
ydata = double(data.image(100:120,100:120,:));
xdata = data.EchoTime';

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Start timing.
tic;

solution = pixelFit(xdata, ydata);

% Stop timing.
toc;

% Plot results.
plot(squeeze(ydata(10,10,:)))
hold
plot(solution(1,10,10)* exp(-xdata'/solution(2,10,10)) + solution(3,10,10))
