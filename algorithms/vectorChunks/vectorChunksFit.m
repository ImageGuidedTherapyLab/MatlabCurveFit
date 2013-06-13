function solution = vectorChunksFit(xdata, ydata)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on', 'MaxIter', 100, 'MaxFunEvals', 100, 'TolFun', 1e-5);
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);

    % Reshape to only one dimension for pixels.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
    solution = ones(2, numberOfPixels);

    % Bounds.
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = 4096 * upperBounds(1,:);
    upperBounds(2,:) = 1500 * upperBounds(2,:);
    lowerBounds = zeros(2, numberOfPixels);
    
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
        solutionChunk{thread} = lsqcurvefit(@vectorT2Decay, solutionChunk{thread}, xdata, ydataChunk{thread}, lowerBoundsChunk{thread}, upperBoundsChunk{thread}, options);
        
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
