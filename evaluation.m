% Reset MATLAB environment.
clear all;
close all;

% Add algorithms.
addpath('algorithms/pixel/');
addpath('algorithms/vector/');

% Load data.
data = load('data/data.mat');
ydata = double(data.image(50:200,30:250,:));
xdata = data.EchoTime';

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Pixel-by-pixel curve fit.
tic();
solutionPixel = pixelFit(xdata, ydata);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata);
processingTimeVector = toc();

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:    % 4.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:  % 4.2f s\n', processingTimeVector);

% Plot results.
figure(1);
plot(squeeze(ydata(10,10,:)));
hold;
plot(solutionPixel(1,10,10)* exp(-xdata'/solutionPixel(2,10,10)) + solutionPixel(3,10,10));
plot(solutionVector(1,10,10)* exp(-xdata'/solutionVector(2,10,10)) + solutionVector(3,10,10));
