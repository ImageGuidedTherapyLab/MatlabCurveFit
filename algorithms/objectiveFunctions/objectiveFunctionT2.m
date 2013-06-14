%% Objective function for lsqcurvefit modeling T2 decay.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier
function [ modelVector, modelJacobian ] = objectiveFunctionT2( solutionParameters, echoTimes )

    %% Initialize.
    % Get number of echos, parameters, and pixels.
    numberOfEchos       = size(echoTimes,1);
    numberOfParameters  = size(solutionParameters,1);
    numberOfPixels      = size(solutionParameters,2);

    % Get each model parameter set.
    amplitudes = solutionParameters(1,:);
    t2Times    = solutionParameters(2,:);

    %% Evaluate objective function.
    % Initalize vector.
    modelVector = zeros(numberOfPixels, numberOfEchos);

    % Build array of model predicted values.
    for echoIndex = 1:numberOfEchos
        modelVector(:,echoIndex) = amplitudes .* exp(-echoTimes(echoIndex) ./ t2Times);
    end

    %% Calculate Jacobian matrix.
    if nargout > 1
        
        % Initialize Jacobian of the function.
        modelJacobian = sparse([], [], [], numberOfPixels*numberOfEchos, numberOfPixels*numberOfParameters, numberOfPixels*numberOfEchos*numberOfParameters);

        for pixelIndex = 1:numberOfPixels
            
            for echoIndex = 1:numberOfEchos
                
                % Index of function components that will be evaluated. The
                % number of function components is the number of pixels
                % times the number of different echo times.
                dataIndex = (echoIndex-1)*numberOfPixels + pixelIndex;
                
                % Index of derivative variables. The number of derivative
                % variables are the number of pixels times the number of
                % fitting parameters. For each function component only one
                % pixel needs to be evaluated. Therefore, the number of
                % calculated derivatives per column equals the number of
                % fitting parameters.
                
                % Parameter 1: Amplitude
                parameterIndex = 1;
                funcIndex = (pixelIndex-1)*numberOfParameters + parameterIndex;
                modelJacobian(dataIndex, funcIndex) = exp(-echoTimes(echoIndex) / t2Times(pixelIndex));

                % Parameter 2: T2
                parameterIndex = 2;
                funcIndex = (pixelIndex-1)*numberOfParameters + parameterIndex;
                modelJacobian(dataIndex, funcIndex) = amplitudes(pixelIndex) * exp(-echoTimes(echoIndex) / t2Times(pixelIndex)) * (echoTimes(echoIndex) / t2Times(pixelIndex)^2);

            end
            
        end
        
    end

end
