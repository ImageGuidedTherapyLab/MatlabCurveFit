%% Parallel pixelwise fitting.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: Chris MacLellan, Florian Maier
function solution = chrisT1Fit(xdata, ydata, solution)

    % Get size of data.
    dataSize = size(ydata);
    numberOfPixels  = dataSize(1)*dataSize(2);
    numberOfSamples = dataSize(3);
        
    % Reshape.
    ydata = reshape(ydata, [numberOfPixels, numberOfSamples]);

    t1sig=@(beta,xx) beta(1)*abs((1-2*exp(-xx/beta(2))));
    
    % Loop over pixels.
    parfor i = 1:numberOfPixels

        warning off;
        
        try
            solution(:,i) = nlinfit( xdata', squeeze(double(ydata(i,:))), t1sig, solution(:,i) );
        catch exception
            fprintf('error: %s', exception.message);
            solution(:,i) = zeros(1,2);
        end
 
    end

    % Reshape solution to format of ydata.
    solution = reshape(solution, [2, dataSize(1), dataSize(2)]);

end
