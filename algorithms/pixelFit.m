%% Parallel pixelwise fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function solution = pixelFit(xdata, ydata, objectiveFunction, solution, bounds)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on' );

    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels  = dataSize(1)*dataSize(2);
    numberOfSamples = dataSize(3);
        
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfSamples]);

    % Bounds.
    lowerBounds(1) = bounds(1,1);
    upperBounds(1) = bounds(1,2);
    lowerBounds(2) = bounds(2,1);
    upperBounds(2) = bounds(2,2);

    % Loop over pixels.
    parfor i = 1:numberOfPixels

        % Fitting.
        solution(:,i) = lsqcurvefit(objectiveFunction, solution(:,i), xdata, ydata(i,:), lowerBounds, upperBounds, options);
 
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
