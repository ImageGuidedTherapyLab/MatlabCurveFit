function [ modelVector, modelJacobian ] = pixelT2Decay( solutionParameters, echoTimes )

    % Get number of echo.
    numberOfEchos = size(echoTimes,1);
    
    % Transpose echo times vector.
    echoTimes = echoTimes';

    % Extract model parameters.
    numberOfParameters = size(solutionParameters,1);
    amplitude = solutionParameters(1);
    decayTime = solutionParameters(2);

    % Build array of model predicted values.
    modelVector = amplitude * exp(-echoTimes / decayTime);
    
    % Calculate Jacobian matrix.
    if nargout > 1
        
        % Jacobian of the function evaluated for solution parameters.
        amplitudePartialDerivatives = exp(-echoTimes ./ decayTime);
        decayPartialDerivatives     = amplitude * exp(-echoTimes ./ decayTime) .* (echoTimes ./ decayTime^2);
        derivativeMatrix = [ amplitudePartialDerivatives decayPartialDerivatives ];
        
        % Create sparse uncouple matrix
        [~, sparseRow] = meshgrid(1:numberOfParameters, 1:numberOfEchos);
        sparseCol = repmat(1:numberOfParameters, numberOfEchos, 1);
        modelJacobian = sparse(sparseRow(:), sparseCol(:), derivativeMatrix(:), numberOfEchos, numberOfParameters);
        
    end
    
end
