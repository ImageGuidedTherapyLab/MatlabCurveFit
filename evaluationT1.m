%% Evaluation script, T1 recovery fitting.
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
data = load('data/dataT1.mat');
ydata = double(data.images);
xdata = double(data.inversionTimes);

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Inital guess for solution.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
initialAmplitude  = max(ydata,[],3);
[~,minimaIndices] = min(ydata,[],3);
initialT1         = xdata(minimaIndices) ./ (1-exp(1));
initialAmplitude  = reshape(initialAmplitude, [1, numberOfPixels]);
initialT1 = reshape(initialT1, [1, numberOfPixels]);
initialGuess = [ initialAmplitude; initialT1 ];

% Define bounds.
bounds = [ 0, 4096; 0, 15000 ];

% Pixel-by-pixel curve fit.
tic();
%solutionPixel = pixelFit(xdata, ydata, @objectiveFunctionT1, initalGuess, bounds);
solutionPixel = reshape(initialGuess, [2, dataSize(1), dataSize(2)]);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, @objectiveFunctionT1, initialGuess, bounds);
processingTimeVector = toc();

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata, @objectiveFunctionT1, initialGuess, bounds);
processingTimeVectorChunks = toc();

fprintf('Processing times:\n');
%fprintf('   pixel-by-pixel, parfor:                     % 4.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                   % 4.2f s\n', processingTimeVector);
fprintf('   vector chunks, piecewise simultaneous fit:  % 4.2f s\n', processingTimeVectorChunks);


% Plot results.
plotRecovery(107,111,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);    % Water
plotRecovery(97,85,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(136,129,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(135,85,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(180,49,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
