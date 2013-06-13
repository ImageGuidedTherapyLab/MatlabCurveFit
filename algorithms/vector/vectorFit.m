% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function solution = vectorFit(xdata, ydata)

    % Set options for fit.
    options = optimset('display', 'iter', 'jacobian', 'on', 'MaxIter', 100, 'MaxFunEvals', 100, 'TolFun', 1e-5);
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);

    % Reshape to avoid for-loops in vectorT2Decay function.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
    solution = ones(2, numberOfPixels);

    % Bounds.
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = 4096 * upperBounds(1,:);
    upperBounds(2,:) = 1500 * upperBounds(2,:);
    lowerBounds = zeros(2, numberOfPixels);
    
    % Fitting.
    solution = lsqcurvefit(@vectorT2Decay, solution, xdata, ydata, lowerBounds, upperBounds, options);

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);
    
end
