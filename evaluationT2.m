%% Evaluation script, T2 decay fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

% Reset MATLAB environment.
clear all;
close all;

% Add paths.
addpath('tools/');
addpath('algorithms/');
addpath('algorithms/objectiveFunctions/');

% Load data.
data = load('data/dataT2.mat');
ydata = double(data.image(50:60,30:40,:));
xdata = data.EchoTime';

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Initial guess for solution.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
initialGuess = ones(2, numberOfPixels);

% Pixel-by-pixel curve fit.
tic();
solutionPixel = pixelFit(xdata, ydata, @objectiveFunctionT2, initialGuess);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, @objectiveFunctionT2, initialGuess);
processingTimeVector = toc();

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata, @objectiveFunctionT2, initialGuess);
processingTimeVectorChunks = toc();

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:                     % 4.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                   % 4.2f s\n', processingTimeVector);
fprintf('   vector chunks, piecewise simultaneous fit:  % 4.2f s\n', processingTimeVectorChunks);


% Plot results.
plotDecay(5,5,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
