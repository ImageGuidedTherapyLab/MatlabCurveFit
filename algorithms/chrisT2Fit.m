%% Parallel pixelwise fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: Chris MacLellan, Florian Maier
function solution = chrisT2Fit(xdata, ydata)

    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels  = dataSize(1)*dataSize(2);
    numberOfSamples = dataSize(3);
        
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfSamples]);

    solution = zeros([2 numberOfPixels]);
    
    % Loop over pixels.
    parfor i = 1:numberOfPixels

        fit = polyfit(xdata', log(squeeze(ydata(i,:))), 1);
        solution(:,i) = [ exp(fit(2)), -1 / fit(1) ];
 
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
