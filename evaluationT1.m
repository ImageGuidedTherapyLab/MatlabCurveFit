% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

% Reset MATLAB environment.
clear all;
close all;

% Add algorithms.
addpath('algorithms/pixel/');
addpath('algorithms/vector/');
addpath('algorithms/vectorChunks/');

% Load data.
data = load('data/dataT1.mat');
ydata = double(data.images(100:150,130:150,:));
xdata = data.inversionTimes;

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Pixel-by-pixel curve fit.
tic();
solutionPixel = pixelFit(xdata, ydata, 'T1');
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, 'T1');
processingTimeVector = toc();

% Vector chunks curve fit.
%tic();
%solutionVectorChunks = vectorChunksFit(xdata, ydata, 'T1');
%processingTimeVectorChunks = toc();

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:                     % 4.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                   % 4.2f s\n', processingTimeVector);
%fprintf('   vector chunks, piecewise simultaneous fit:  % 4.2f s\n', processingTimeVectorChunks);


% Plot results.
%plotDecay(60,80,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks,'T1');
