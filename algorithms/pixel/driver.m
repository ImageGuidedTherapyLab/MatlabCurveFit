clear all
close all
load ../../data/data.mat

options = optimset('display','off','jacobian','on','MaxIter',40,'MaxFunEvals',40,'TolFun',1e-6);
ydata = double(image(50:200,30:250,:));
xdata = EchoTime'

matlabpool open;


tic

% Get size of data.
dataSize = size(ydata);
numberOfPixels = dataSize(1)*dataSize(2);
numberOfEchos  = dataSize(3);

% Reshape to avoid nested for-loops.
ydata = reshape(ydata, [numberOfPixels, numberOfEchos]);
solution = zeros(3, numberOfPixels);

lowerBounds(1) = 0;
lowerBounds(2) = 0;
lowerBounds(3) = 0;

upperBounds(1) = 200;
upperBounds(2) =  60;
upperBounds(3) = 100;

% Loop over pixels.
parfor i = 1:numberOfPixels
    
    echoTimes = xdata;
    echoData  = ydata(i,:);
    
    x0 = rand(3,1);

    
    [x,resnorm,residual,exitflag, output] = lsqcurvefit( @t2decay, x0, echoTimes, echoData, lowerBounds, upperBounds, options);
    solution(:,i) = x(:);
    
end

% Reshape data to previous format.
ydata = reshape(ydata, dataSize);
solution = reshape(solution, [3, dataSize(1), dataSize(2)]);

toc

plot(squeeze(ydata(60,80,:)))
hold
plot(solution(1,60,80)* exp(-EchoTime/solution(2,60,80)) + solution(3,60,80))

matlabpool close;
