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
ydata = double(data.image(100:179,50:129,:));
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
solutionPixel = pixelFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
processingTimePixel = toc();

% Vector curve fit.
tic();
solutionVector = vectorFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
processingTimeVector = toc();

% Vector chunks curve fit.
tic();
solutionVectorChunks = vectorChunksFit(xdata, ydata, @objectiveFunctionT2, initialGuess, bounds);
processingTimeVectorChunks = toc();

% Chris' curve fit.
tic();
solutionChris = chrisT2Fit(xdata, ydata);
processingTimeChris = toc();

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
fprintf('   Chris'' way:                                % 7.2f s\n', processingTimeChris);

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

echoTimes  = repmat( reshape( xdata, [1 1 12]), [size(ydata,1) size(ydata,2) 1]);
amplitudes = repmat( reshape( solutionPixel(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t2times    = repmat( reshape( solutionPixel(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f\n', error);

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

amplitudes = repmat( reshape( solutionVector(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t2times    = repmat( reshape( solutionVector(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f\n', error);

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

amplitudes = repmat( reshape( solutionVectorChunks(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t2times    = repmat( reshape( solutionVectorChunks(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f\n', error);

figure;
imagesc(squeeze(solutionChris(2,:,:)));
title('T_2 map (Chris'' method)');
colormap jet;
xlabel('x / px');
ylabel('y / px');
caxis([0 100]);
h = colorbar;
set(get(h,'ylabel'),'String', 'T2 / ms');
snapnow;

amplitudes = repmat( reshape( solutionChris(1,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);
t2times    = repmat( reshape( solutionChris(2,:,:), [size(ydata,1) size(ydata,2) 1]), [1 1 size(xdata,1)]);

valid = ((~isnan(amplitudes)) & (~isnan(t2times)) & (amplitudes > 0) & (t2times > 0));
error = rmse( ydata(valid), amplitudes(valid) .* exp( -echoTimes(valid) ./ t2times(valid)) );
fprintf('RMSE = %.3f\n', error);

% Plot results.
plotDecay( 26, 61, xdata, ydata, solutionPixel, solutionVector, solutionVectorChunks, solutionChris ); % Liver
plotDecay(  5, 78, xdata, ydata, solutionPixel, solutionVector, solutionVectorChunks, solutionChris ); % Fat
plotDecay( 80, 45, xdata, ydata, solutionPixel, solutionVector, solutionVectorChunks, solutionChris ); % Aorta
plotDecay( 73, 59, xdata, ydata, solutionPixel, solutionVector, solutionVectorChunks, solutionChris ); % Bone
