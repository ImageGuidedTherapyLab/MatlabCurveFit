%% Plot T1 recovery.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function plotRecovery(x,y,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks)

    % Plot results.
    pixel = [y,x];
    figure;
    plot(xdata,squeeze(ydata(pixel(1),pixel(2),:)),'+',...
         1:max(xdata),solutionPixel(1,pixel(1),pixel(2))* abs(1 - 2 * exp(-[1:max(xdata)]/solutionPixel(2,pixel(1),pixel(2)))),'g-',...
         1:max(xdata),solutionVector(1,pixel(1),pixel(2))* abs(1 - 2 * exp(-[1:max(xdata)]/solutionVector(2,pixel(1),pixel(2)))),'b-',...
         1:max(xdata),solutionVectorChunks(1,pixel(1),pixel(2))* abs(1 - 2 * exp(-[1:max(xdata)]/solutionVectorChunks(2,pixel(1),pixel(2)))),'r-');

end
