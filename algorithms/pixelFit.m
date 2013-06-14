%% Parallel pixelwise fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function solution = pixelFit(xdata, ydata, objectiveFunction, solution)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on' );

    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);
        
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);

    % Bounds for amplitude
    lowerBounds(1) = 0;
    upperBounds(1) = 4096;

    % Bounds for T1 time
    lowerBounds(2) = 0;
    upperBounds(2) = 5000;

    % Loop over pixels.
    parfor i = 1:numberOfPixels

        echoTimes = xdata;
        echoData  = ydata(i,:);

        x0 = solution(:,i);
        
        % Fitting.
        solution(:,i) = lsqcurvefit(objectiveFunction, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
 
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
