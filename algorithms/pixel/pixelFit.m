function solution = pixelFit(xdata, ydata)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on', 'MaxIter', 40, 'MaxFunEvals', 40, 'TolFun', 1e-3 );

    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);

    % Reshape to avoid nested for-loops.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
    solution = ones(2, numberOfPixels);

    % Bounds for amplitude
    lowerBounds(1) = 0;
    upperBounds(1) = 4096;

    % Bounds for T2 time
    lowerBounds(2) = 0;
    upperBounds(2) = 1500;

    % Loop over pixels.
    parfor i = 1:numberOfPixels

        echoTimes = xdata;
        echoData  = ydata(i,:);

        x0 = solution(:,i);
        solution(:,i) = lsqcurvefit( @pixelT2Decay, x0, echoTimes, echoData, lowerBounds, upperBounds, options);

    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
