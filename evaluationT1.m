%% Evaluation script, T1 recovery fitting.

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
initialT1         = xdata(minimaIndices);
initialT1         = initialT1 ./ (1-1/exp(1));
initialAmplitude  = reshape(initialAmplitude, [1, numberOfPixels]);
initialT1 = reshape(initialT1, [1, numberOfPixels]);
initialGuess = [ initialAmplitude; initialT1 ];

% Define bounds.
bounds = [ 0, 4096; 0, 15000 ];

% Pixel-by-pixel curve fit.
tic();
%solutionPixel = pixelFit(xdata, ydata, @objectiveFunctionT1, initialGuess, bounds);
%solutionPixel = reshape(initialGuess, [2, dataSize(1), dataSize(2)]);
solutionPixel = chrisT1Fit(xdata, ydata, initialGuess);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, @objectiveFunctionT1, initialGuess, bounds);
processingTimeVector = toc();

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata, @objectiveFunctionT1, initialGuess, bounds);
processingTimeVectorChunks = toc();

figure;
imagesc(ydata(:,:,1));
title('magnitude image');
colormap gray;
xlabel('x / px');
ylabel('y / px');
h = colorbar;
set(get(h,'ylabel'),'String', 'signal magnitude / a.u.');
snapnow;

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:                    % 7.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                  % 7.2f s\n', processingTimeVector);
fprintf('   vector chunks, piecewise simultaneous fit: % 7.2f s\n', processingTimeVectorChunks);

figure;
imagesc(squeeze(solutionPixel(2,:,:)));
title('T_1 map (pixel-by-pixel, parfor)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 2000]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T1 / ms');
snapnow;

inversionTimes = repmat( reshape( xdata, [1 1 size(xdata,1)]), [size(ydata,1) size(ydata,2) 1]);
amplitudes     = repmat( reshape( solutionPixel(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t1times        = repmat( reshape( solutionPixel(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t1times)) & (amplitudes > 0) & (t1times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -inversionTimes(valid) ./ t1times(valid)) );
fprintf('RMSE = %.3f', error);

figure;
imagesc(squeeze(solutionVector(2,:,:)));
title('T_1 map (vector, simultaneous fit)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 2000]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T1 / ms');
snapnow;

amplitudes = repmat( reshape( solutionVector(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t1times    = repmat( reshape( solutionVector(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t1times)) & (amplitudes > 0) & (t1times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -inversionTimes(valid) ./ t1times(valid)) );
fprintf('RMSE = %.3f', error);

figure;
imagesc(squeeze(solutionVectorChunks(2,:,:)));
title('T_1 map (vector chunks, piecewise simultaneous fit)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 2000]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T1 / ms');
snapnow;

amplitudes = repmat( reshape( solutionVectorChunks(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t1times    = repmat( reshape( solutionVectorChunks(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t1times)) & (amplitudes > 0) & (t1times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -inversionTimes(valid) ./ t1times(valid)) );
fprintf('RMSE = %.3f', error);

% Plot results.
plotRecovery(107,111,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);    % Water
plotRecovery(97,85,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(136,129,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(135,85,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
plotRecovery(180,49,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);
