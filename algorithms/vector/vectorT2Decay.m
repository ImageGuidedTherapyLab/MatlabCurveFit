% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function [ModelVector,ModelJacobian] = vectorT2Decay(Solution,xdata)

    % get number of echo
    NumEcho = size(xdata,1);

    % objective function values at Solution
    ModelVector = zeros(size(Solution,2),NumEcho);

    % extract model parameters individually
    Amplitude = Solution(1,:);
    Decay     = Solution(2,:);
    Nparam    = 2; 

    % build array of model predicted values
    for iii = 1:NumEcho
        echotime = xdata(iii,1);
        ModelVector(:,iii) =  Amplitude .* exp(-echotime./Decay);
    end
    
    if nargout > 1   % two output arguments
        
        NumberRowsPerEcho = size(Solution,2);
        NumberRowsTotal   = size(Solution,2) * size(xdata,1);
   
        % allocate dense matrix for derivative
        DerivativeMatrix = zeros(NumberRowsTotal,Nparam);
   
        % Jacobian of the function evaluated at Solution
        % For jacobian of matrix function and variables, see:
        % http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
        for iii = 1:size(xdata,1)
            echotime = xdata(iii,1);
            AmplitudePartialDerivative = exp(-echotime./Decay(:));
            DecayPartialDerivative     = Amplitude(:).*exp(-echotime./Decay(:)).*(echotime./Decay(:).^2);
            DerivativeMatrix( ((iii-1)*NumberRowsPerEcho+1):iii*NumberRowsPerEcho,:) = [AmplitudePartialDerivative, DecayPartialDerivative ];
        end
        
        % create sparse uncouple matrix
        [~,SparseRow] = meshgrid([1:Nparam],[1:NumberRowsTotal]);
        TmpCol = reshape([1:Nparam*NumberRowsPerEcho],Nparam,NumberRowsPerEcho)';
        SparseCol = repmat(TmpCol,NumEcho,1);
        ModelJacobian = sparse(SparseRow(:),SparseCol(:),DerivativeMatrix(:),NumberRowsTotal,Nparam*NumberRowsPerEcho);
        
    end

end
