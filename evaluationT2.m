%% Evaluation script, T2 decay fitting.

% Reset MATLAB environment.
clear all;
close all;

% Add paths.
addpath('tools/');
addpath('algorithms/');
addpath('algorithms/objectiveFunctions/');

% Load data.
data = load('data/dataT2.mat');
ydata = double(data.image);
xdata = data.EchoTime';

% Create pool for parallel processing.
if matlabpool('size') == 0
    matlabpool('open');
end

% Initial guess for solution.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
initialGuess = ones(2, numberOfPixels);

% Define bounds.
bounds = [ 0, 4096; 0, 2500 ];

% Pixel-by-pixel curve fit.
tic();
%solutionPixel = pixelFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
%solutionPixel = reshape(initialGuess, [2, dataSize(1), dataSize(2)]);
solutionPixel = chrisT2Fit(xdata, ydata);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
processingTimeVector = toc();

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
processingTimeVectorChunks = toc();

figure;
imagesc(ydata(:,:,1));
title('magnitude image');
colormap gray;
xlabel('x / px');
ylabel('y / px');
caxis([0 800]);
h = colorbar;
set(get(h,'ylabel'),'String', 'signal magnitude / a.u.');
snapnow;

fprintf('Processing times:\n');
fprintf('   pixel-by-pixel, parfor:                    % 7.2f s\n', processingTimePixel);
fprintf('   vector, simultaneous fit:                  % 7.2f s\n', processingTimeVector);
fprintf('   vector chunks, piecewise simultaneous fit: % 7.2f s\n', processingTimeVectorChunks);

figure;
imagesc(squeeze(solutionPixel(2,:,:)));
title('T_2 map (pixel-by-pixel, parfor)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 100]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T2 / ms');
snapnow;

echoTimes  = repmat( reshape( xdata, [1 1 12]), [256 256 1]);
amplitudes = repmat( reshape( solutionPixel(1,:,:), [256 256 1]), [1 1 12]);
t2times    = repmat( reshape( solutionPixel(2,:,:), [256 256 1]), [1 1 12]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f', error);

figure;
imagesc(squeeze(solutionVector(2,:,:)));
title('T_2 map (vector, simultaneous fit)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 100]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T2 / ms');
snapnow;

amplitudes = repmat( reshape( solutionVector(1,:,:), [256 256 1]), [1 1 12]);
t2times    = repmat( reshape( solutionVector(2,:,:), [256 256 1]), [1 1 12]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f', error);

figure;
imagesc(squeeze(solutionVectorChunks(2,:,:)));
title('T_2 map (vector chunks, piecewise simultaneous fit)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 100]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T2 / ms');
snapnow;

amplitudes = repmat( reshape( solutionVectorChunks(1,:,:), [256 256 1]), [1 1 12]);
t2times    = repmat( reshape( solutionVectorChunks(2,:,:), [256 256 1]), [1 1 12]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f', error);

% Plot results.
plotDecay(90,128,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);    % Liver
plotDecay(40,103,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);    % Fat
plotDecay(134,143,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);   % Aorta
plotDecay(231,132,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks);   % Muscle
