%% Full vector fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function solution = vectorFit(xdata, ydata, objectiveFunction, solution)

    % Set options for fit.
    options = optimset('display', 'iter', 'jacobian', 'on');
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);
    
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);

    % Bounds.
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = 4096 * upperBounds(1,:);
    upperBounds(2,:) = 5000 * upperBounds(2,:);
    lowerBounds = zeros(2, numberOfPixels);
    
    % Fitting.
    solution = lsqcurvefit(objectiveFunction, solution, xdata, ydata, lowerBounds, upperBounds, options);

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);
    
end
