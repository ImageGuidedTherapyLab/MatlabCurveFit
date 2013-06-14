% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function [ modelVector, modelJacobian ] = pixelT1Recovery( solutionParameters, inversionTimes )

    % Get number of inversion times.
    numberOfInversionTimes = size(inversionTimes,1);
    
    % Reshape inversion times vector.
    inversionTimes = reshape(inversionTimes, [1,numberOfInversionTimes]);

    % Extract model parameters.
    numberOfParameters = size(solutionParameters,1);
    amplitude = solutionParameters(1);
    t1 = solutionParameters(2);

    % Build array of model predicted values.
    modelVector = amplitude * abs(1 - 2*exp(-inversionTimes / t1));
    
    % Calculate Jacobian matrix.
    if nargout > 1
        
        % Jacobian of the function evaluated for solution parameters.
        amplitudePartialDerivatives = abs(1 - 2*exp(-inversionTimes / t1));
        t1PartialDerivatives        = sign(1 - 2*exp(-inversionTimes ./ t1)) .* amplitude .* (1 - 2*exp(-inversionTimes / t1)) .* (-2*exp(-inversionTimes / t1)) .* (inversionTimes ./ t1^2);
        derivativeMatrix = [ amplitudePartialDerivatives t1PartialDerivatives ];
        
        % Create sparse uncouple matrix
        [~, sparseRow] = meshgrid(1:numberOfParameters, 1:numberOfInversionTimes);
        sparseCol = repmat(1:numberOfParameters, numberOfInversionTimes, 1);
        modelJacobian = sparse(sparseRow(:), sparseCol(:), derivativeMatrix(:), numberOfInversionTimes, numberOfParameters);

    end
    
end
