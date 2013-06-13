function [ModelVector,ModelJacobian] = t2decay(Solution,xdata)

    % get number of echo
    NumEcho = size(xdata,1);

    % objective function values at Solution
    ModelVector = zeros(1,NumEcho);

    % extract model parameters individually
    Amplitude = Solution(1);
    Decay     = Solution(2);
    Offset    = Solution(3);
    Nparam    = 3; 

    % build array of model predicted values
    for iii = 1:NumEcho
        echotime = xdata(iii);
        ModelVector(iii) =  Amplitude * exp(-echotime / Decay) + Offset;
    end
    
    if nargout > 1   % two output arguments
        NumberRowsPerEcho = 1;
        NumberRowsTotal   = NumEcho;
        
        % allocate dense matrix for derivative
        DerivativeMatrix = zeros(NumberRowsTotal,3);
   
        % Jacobian of the function evaluated at Solution
        % For jacobian of matrix function and variables, see:
        % http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
        for iii = 1:size(xdata,1)
            echotime = xdata(iii);
            AmplitudePartialDerivative = exp(-echotime/Decay);
            DecayPartialDerivative     = Amplitude * exp(-echotime/Decay) * (echotime/Decay^2);
            OffsetPartialDerivative    = 1;
            DerivativeMatrix(iii,:) = [AmplitudePartialDerivative, DecayPartialDerivative, OffsetPartialDerivative ];
        end
        
        % Create sparse uncouple matrix
        [NotUsed,SparseRow] = meshgrid([1:Nparam],[1:NumberRowsTotal]);
        TmpCol = 1:Nparam;
        SparseCol = repmat(TmpCol,NumEcho,1);
        ModelJacobian = sparse(SparseRow(:),SparseCol(:),DerivativeMatrix(:),NumberRowsTotal,Nparam);
    end
    
end