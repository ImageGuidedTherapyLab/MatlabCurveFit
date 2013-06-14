% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function solution = vectorChunksFit(xdata, ydata, mode)

    % Set options for fit.
    options = optimset('display', 'off', 'jacobian', 'on');
    
    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels = dataSize(1)*dataSize(2);
    numberOfEchos  = dataSize(3);

    switch (mode)
        case 'T1'
            % Guess inital parameters.
            initialAmplitude  = max(ydata,[],3);
            [~,minimaIndices] = min(ydata,[],3);
            initialT1         = xdata(minimaIndices);
            initialAmplitude = reshape(initialAmplitude, [1, numberOfPixels]);
            initialT1 = reshape(initialT1, [1, numberOfPixels]);
            solution = [ initialAmplitude; initialT1 ];
        case 'T2'
            solution = ones(2, numberOfPixels);
    end
            
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);

    % Bounds.
    upperBounds = ones(2, numberOfPixels);
    upperBounds(1,:) = 4096 * upperBounds(1,:);
    upperBounds(2,:) = 5000 * upperBounds(2,:);
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
         switch (mode)
             case 'T1'
                solutionChunk{thread} = lsqcurvefit(@objectiveFunctionT1, solutionChunk{thread}, xdata, ydataChunk{thread}, lowerBoundsChunk{thread}, upperBoundsChunk{thread}, options);
             case 'T2'
                solutionChunk{thread} = lsqcurvefit(@objectiveFunctionT2, solutionChunk{thread}, xdata, ydataChunk{thread}, lowerBoundsChunk{thread}, upperBoundsChunk{thread}, options);
         end
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
