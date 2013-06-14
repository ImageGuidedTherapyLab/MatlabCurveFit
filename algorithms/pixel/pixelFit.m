% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function solution = pixelFit(xdata, ydata, mode)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on' );

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
        
        switch(mode)
            case 'T1'
                solution(:,i) = lsqcurvefit( @pixelT1Recovery, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
            case 'T2'
                solution(:,i) = lsqcurvefit( @pixelT2Decay, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
        end

    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
