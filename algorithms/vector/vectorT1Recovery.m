% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function [ modelVector, modelJacobian ] = vectorT1Recovery( solutionParameters, inversionTimes )

    %% Initialize.
    % Get number of inversion times, parameters, and pixels.
    numberOfInversionTimes = size(inversionTimes,1);
    numberOfParameters     = size(solutionParameters,1);
    numberOfPixels         = size(solutionParameters,2);

    % Get each model parameter set.
    amplitudes = solutionParameters(1,:);
    t1Times    = solutionParameters(2,:);

    %% Evaluate objective function.
    % Initalize vector.
    modelVector = zeros(numberOfPixels, numberOfInversionTimes);

    % Build array of model predicted values.
    for inversionIndex = 1:numberOfInversionTimes
        modelVector(:,inversionIndex) = amplitudes .* abs(1 - 2.*exp(-inversionTimes(inversionIndex) ./ t1Times));
    end

    %% Calculate Jacobian matrix.
    if nargout > 1
        
        % Initialize Jacobian of the function.
        modelJacobian = sparse([], [], [], numberOfPixels*numberOfInversionTimes, numberOfPixels*numberOfParameters, numberOfPixels*numberOfInversionTimes*numberOfParameters);

        for pixelIndex = 1:numberOfPixels
            
            for inversionIndex = 1:numberOfInversionTimes
                
                % Index of function components that will be evaluated. The
                % number of function components is the number of pixels
                % times the number of different inversion times.
                dataIndex = (inversionIndex-1)*numberOfPixels + pixelIndex;
                
                % Index of derivative variables. The number of derivative
                % variables are the number of pixels times the number of
                % fitting parameters. For each function component only one
                % pixel needs to be evaluated. Therefore, the number of
                % calculated derivatives per column equals the number of
                % fitting parameters.
                
                % Parameter 1: Amplitude
                parameterIndex = 1;
                funcIndex = (pixelIndex-1)*numberOfParameters + parameterIndex;
                modelJacobian(dataIndex, funcIndex) = abs(1 - 2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex)));
                
                % Parameter 2: T1
                parameterIndex = 2;
                funcIndex = (pixelIndex-1)*numberOfParameters + parameterIndex;
                modelJacobian(dataIndex, funcIndex) = sign(1 - 2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex))) * amplitudes(pixelIndex) * (1 - 2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex))) * (-2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex))) * (inversionTimes(inversionIndex) / t1Times(pixelIndex)^2);
                
            end
            
        end
        
    end

end
