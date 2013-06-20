%% Objective function for lsqcurvefit modeling T1 recovery.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function [ modelVector, modelJacobian ] = objectiveFunctionT1( solutionParameters, inversionTimes )

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
    for pixelIndex = 1:numberOfPixels
        
        for inversionIndex = 1:numberOfInversionTimes
            
            modelVector(pixelIndex,inversionIndex) = amplitudes(pixelIndex) * abs(1 - 2 * exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex)));
        
        end
        
    end

    %% Calculate Jacobian matrix.
    if nargout > 1
        
        % Initialize vectors for Jacobian of the function.
        jacobianRowIndices = zeros(1,numberOfPixels*numberOfInversionTimes*numberOfParameters);
        jacobianColIndices = zeros(1,numberOfPixels*numberOfInversionTimes*numberOfParameters);
        jacobianEntries    = zeros(1,numberOfPixels*numberOfInversionTimes*numberOfParameters);
        i = 1;
        
        for pixelIndex = 1:numberOfPixels
            
            for inversionIndex = 1:numberOfInversionTimes
                
                % Parameter 1: Amplitude
                parameterIndex = 1;
                jacobianRowIndices(i)   = (inversionIndex-1)*numberOfPixels + pixelIndex;
                jacobianColIndices(i)   = (pixelIndex-1)*numberOfParameters + parameterIndex;
                jacobianEntries(i)      = abs(1 - 2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex)));
                
                % Parameter 2: T1
                parameterIndex = 2;
                jacobianRowIndices(i+1) = (inversionIndex-1)*numberOfPixels + pixelIndex;
                jacobianColIndices(i+1) = (pixelIndex-1)*numberOfParameters + parameterIndex;
                jacobianEntries(i+1)    = sign(1 - 2*exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex))) * amplitudes(pixelIndex) * -2 * exp(-inversionTimes(inversionIndex) / t1Times(pixelIndex)) * (inversionTimes(inversionIndex) / t1Times(pixelIndex)^2);
                
                % Update Jacobian indices index.
                i = i + 2;
                
            end
            
        end

        % Create sparse matrix.
        modelJacobian = sparse(jacobianRowIndices, jacobianColIndices, jacobianEntries);
        
    end

end
