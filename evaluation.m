% Reset MATLAB environment.
clear all;
close all;

% Add algorithms.
addpath('algorithms/pixel/');
addpath('algorithms/vector/');
addpath('algorithms/vectorChunks/');

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

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata);
processingTimeVectorChunks = toc();

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:                     % 4.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                   % 4.2f s\n', processingTimeVector);
fprintf('   vector chunks, piecewise simultaneous fit:  % 4.2f s\n', processingTimeVectorChunks);


% Plot results.
plotDecay(60,80,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
