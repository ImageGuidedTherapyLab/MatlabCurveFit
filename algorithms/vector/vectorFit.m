function solution = vectorFit(xdata, ydata)

    % Set options for fit.
    options = optimset('display', 'iter', 'jacobian', 'on', 'MaxIter', 100, 'MaxFunEvals', 100, 'TolFun', 1e-6);
    
    % Get size of data.
    dataSize = size(ydata);

    % Allocate space
    % Array shape wrt jacobian of matrix function and variables, see:
    % http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
    x0 = ones(2, dataSize(1), dataSize(2));
    
    % Bounds.
    upperBounds = ones(2, dataSize(1), dataSize(2));
    upperBounds(1,:,:) = 4096 * upperBounds(1,:,:);
    upperBounds(2,:,:) = 1500 * upperBounds(2,:,:);    
    lowerBounds = zeros(2, dataSize(1), dataSize(2));
    
    % Fitting.
    solution = lsqcurvefit(@vectorT2Decay, x0, xdata, ydata, lowerBounds, upperBounds, options);

end
