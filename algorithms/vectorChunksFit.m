%% Parallel piecewise vector fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function solution = vectorChunksFit(xdata, ydata, objectiveFunction, solution, bounds)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on');
    
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
    
    % Get number of threads.
    numberOfThreads = matlabpool('size');
    
    % Process one chunk per thread.
    parfor thread = 1:numberOfThreads
        
        % Get bounds of chunk.
        chunkSize = ceil(numberOfPixels / numberOfThreads);
        firstPixel = 1 + (thread-1) * chunkSize;
        lastPixel  = min([firstPixel + chunkSize - 1, numberOfPixels]);

        % Get chunks.
        ydataChunk{thread} = ydata(firstPixel:lastPixel,:);
        lowerBoundsChunk{thread} = lowerBounds(:,firstPixel:lastPixel);
        upperBoundsChunk{thread} = upperBounds(:,firstPixel:lastPixel);
        solutionChunk{thread} = solution(:,firstPixel:lastPixel);
        
        % Fitting.
        solutionChunk{thread} = lsqcurvefit(objectiveFunction, solutionChunk{thread}, xdata, ydataChunk{thread}, lowerBoundsChunk{thread}, upperBoundsChunk{thread}, options);

    end
    
    % Combine results.
    for thread = 1:numberOfThreads
        
        % Get bounds of chunk.
        chunkSize = ceil(numberOfPixels / numberOfThreads);
        firstPixel = 1 + (thread-1) * chunkSize;
        lastPixel  = min([firstPixel + chunkSize - 1, numberOfPixels]);
    
        % Copy data back to full solution.
        solution(:,firstPixel:lastPixel) = solutionChunk{thread};
        
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);
    
end
