% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function solution = vectorFit(xdata, ydata, mode)

    % Set options for fit.
    options = optimset('display', 'iter', 'jacobian', 'on', 'MaxIter', 100, 'MaxFunEvals', 100, 'TolFun', 1e-7);
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);
    
    % Guess inital parameters.
    initialAmplitude  = max(ydata,[],3);
    [~,minimaIndices] = min(ydata,[],3);
    initialT1         = xdata(minimaIndices);

    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
    initialAmplitude = reshape(initialAmplitude, [1, numberOfPixels]);
    initialT1 = reshape(initialT1, [1, numberOfPixels]);
    solution = [ initialAmplitude; initialT1 ];

    % Bounds.
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = 4096 * upperBounds(1,:);
    upperBounds(2,:) = 5000 * upperBounds(2,:);
    lowerBounds = zeros(2, numberOfPixels);
    
    % Fitting.
    switch (mode)
        case 'T1'
            solution = lsqcurvefit(@vectorT1Recovery, solution, xdata, ydata, lowerBounds, upperBounds, options);
        case 'T2'
            solution = lsqcurvefit(@vectorT2Decay, solution, xdata, ydata, lowerBounds, upperBounds, options);
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);
    
end
