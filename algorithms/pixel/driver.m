clear all
close all
load ../../data/data.mat

options = optimset('display','iter','jacobian','on','MaxIter',40,'MaxFunEvals',40,'TolFun',1e-6);
ydata = double(image(50:200,30:250,:));
xdata = EchoTime'

tic

% Get size of data.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
numberOfEchos  = dataSize(3);

% Reshape to avoid nested for-loops.
ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);

% Loop over pixels.
for i = 1:numberOfPixels
    
    echoTimes = xdata;
    echoData  = ydata(i,:);
    
    x0 = rand(3,1);

    lowerBounds(1) = 0;
    lowerBounds(2) = 0;
    lowerBounds(3) = 0;
    
    upperBounds(1) = 200;
    upperBounds(2) =  60;
    upperBounds(3) = 100;
    
    [x,resnorm,residual,exitflag, output] = lsqcurvefit( @t2decay, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
    
end

% Reshape data to previous format.
ydata = reshape(ydata, dataSize);

toc

plot(squeeze(ydata(60,80,:)))
hold
plot(x(1,60,80)* exp(-EchoTime./x(2,60,80)) + x(3,60,80))
