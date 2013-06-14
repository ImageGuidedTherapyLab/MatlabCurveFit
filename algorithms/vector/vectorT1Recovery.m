% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function [ modelVector, modelJacobian ] = vectorT1Recovery( solutionParameters, inversionTimes )

    % Get number of inversion times.
    numberOfInversionTimes = size(inversionTimes,1);

    % Reshape inversion times vector.
    inversionTimes = reshape(inversionTimes, [1,numberOfInversionTimes]);

    % Get number of pixels.
    numberOfPixels = size(solutionParameters,2);
    
    % Objective function values at Solution
    modelVector = zeros(numberOfPixels, numberOfInversionTimes);

    % Extract model parameters.
    numberOfParameters = size(solutionParameters,1);
    amplitudes = solutionParameters(1,:);
    t1Times    = solutionParameters(2,:);

    % Build array of model predicted values
    for i = 1:numberOfInversionTimes
        inversionTime = inversionTimes(i);
        modelVector(:,i) = amplitudes .* abs(1 - 2.*exp(-inversionTime ./ t1Times));
    end

    % Calculate Jacobian matrix.
    if nargout > 1
        
        numberOfRows = numberOfPixels * numberOfInversionTimes;
   
        % Allocate dense matrix for derivative
        derivativeMatrix = zeros(numberOfRows, numberOfParameters);
   
        % Jacobian of the function evaluated at Solution
        % For jacobian of matrix function and variables, see:
        % http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
        for i = 1:numberOfInversionTimes
            inversionTime = inversionTimes(i);
            amplitudePartialDerivative = abs(1 - 2.*exp(-inversionTime ./ t1Times));
            t1PartialDerivative        = sign(1 - 2.*exp(-inversionTime ./ t1Times)) .* amplitudes .* (1 - 2.*exp(-inversionTime ./ t1Times)) .* (-2.*exp(-inversionTime ./ t1Times)) .* (inversionTime ./ t1Times.^2);
            derivativeMatrix( ((i-1)*numberOfPixels+1):i*numberOfPixels,:) = [ amplitudePartialDerivative', t1PartialDerivative' ];
        end
        
        % create sparse uncouple matrix
        [~,sparseRow] = meshgrid(1:numberOfParameters, 1:numberOfRows);
        tmpCol = reshape(1:(numberOfParameters*numberOfPixels),numberOfPixels,numberOfParameters);
        sparseCol = repmat(tmpCol,numberOfInversionTimes,1);
        modelJacobian = sparse(sparseRow(:),sparseCol(:),derivativeMatrix(:),numberOfRows,numberOfParameters*numberOfPixels);
        
    end

end
