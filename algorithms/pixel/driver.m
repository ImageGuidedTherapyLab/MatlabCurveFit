% Reset MATLAB environment.
clear all
close all

% Start timing.
tic;

% Load data.
load('../../data/data.mat');
ydata = double(image(100:120,100:120,:));
xdata = EchoTime';

% Set options for fit.
options = optimset('display', 'off', 'jacobian', 'on', 'MaxIter', 40, 'MaxFunEvals', 40, 'TolFun', 1e-6 );

% Create pool for parallel processing.
matlabpool open;

% Get size of data.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
numberOfEchos  = dataSize(3);

% Reshape to avoid nested for-loops.
ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
solution = rand(3, numberOfPixels);

% Bounds for amplitude
lowerBounds(1) = 0;
upperBounds(1) = 200;

% Bounds for T2 time
lowerBounds(2) = 0;
upperBounds(2) = 60;

% Bounds for offset
lowerBounds(3) = 0;
upperBounds(3) = 100;

% Loop over pixels.
parfor i = 1:numberOfPixels
    
    echoTimes = xdata;
    echoData  = ydata(i,:);
    
    x0 = solution(:,i);

    [x,resnorm,residual,exitflag, output] = lsqcurvefit( @t2decay, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
    
    solution(:,i) = x(:);
    
end

% Reshape data to previous format.
ydata = reshape(ydata, dataSize);
solution = reshape(solution, [3, dataSize(1), dataSize(2)]);

% Close MATLAB pool
matlabpool close;

% Stop timing.
toc;

% Plot results.
plot(squeeze(ydata(10,10,:)))
hold
plot(solution(1,10,10)* exp(-EchoTime/solution(2,10,10)) + solution(3,10,10))
