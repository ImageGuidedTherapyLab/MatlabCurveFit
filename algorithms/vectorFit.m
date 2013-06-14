%% Full vector fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function solution = vectorFit(xdata, ydata, objectiveFunction, solution, bounds)

    % Set options for fit.
    options = optimset('display', 'iter', 'jacobian', 'on');
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels  = dataSize(1)*dataSize(2);
    numberOfSamples = dataSize(3);
    
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfSamples]);

    % Bounds.
    lowerBounds = ones(2, numberOfPixels);
    lowerBounds(1,:) = bounds(1,1) * lowerBounds(1,:);
    lowerBounds(2,:) = bounds(2,1) * lowerBounds(2,:);
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = bounds(1,2) * upperBounds(1,:);
    upperBounds(2,:) = bounds(2,2) * upperBounds(2,:);
    
    % Fitting.
    solution = lsqcurvefit(objectiveFunction, solution, xdata, ydata, lowerBounds, upperBounds, options);

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);
    
end
