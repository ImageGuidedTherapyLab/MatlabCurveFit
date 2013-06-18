%% Plot T1 recovery.
% Copyright (c) The University of Texas MD Anderson Cancer Center, 2013
% Authors: David Fuentes, Florian Maier

function plotRecovery(x,y,xdata,ydata,solutionPixel,solutionVector,solutionVectorChunks)

    % Init.
    pixel    = [ y, x ];
    measured = squeeze( ydata( pixel(1), pixel(2), : ));
    t        = 0:1.1*max(xdata);
    func     = @(t,param) param(1) * abs(1 - 2 * exp(-t ./ param(2)));
    
    % Plot results.
    figure;
    imagesc(ydata(:,:,1));
    title('magnitude image');
    colormap gray;
    xlabel('x / px');
    ylabel('y / px');
    h = colorbar;
    set(get(h,'ylabel'),'String', 'signal magnitude / a.u.');
    hold on;
    plot( x, y, 'y+', 'MarkerSize', 10);
    snapnow;
    
    figure;
    plot(xdata ,measured, 'k+',...
         t, func( t, solutionPixel(:,pixel(1),pixel(2)) ), 'r-', ...
         t, func( t, solutionVector(:,pixel(1),pixel(2)) ),'g-', ...
         t, func( t, solutionVectorChunks(:,pixel(1),pixel(2)) ),'b--')
    title(sprintf('T_1 recovery at pixel (%d/%d)', x, y));
    xlabel('echo time / ms');
    ylabel('signal / a.u.');
    legend('measurement', 'pixel-by-pixel fit', 'simultaneous fit', 'piecewise simultaneous fit', 'Location', 'SouthEast');
    snapnow;
    
    % RMSE calculations.
    rmsePixel        = rmse( measured, func( xdata, solutionPixel(:,pixel(1),pixel(2))));
    rmseVector       = rmse( measured, func( xdata, solutionVector(:,pixel(1),pixel(2))));
    rmseVectorChunks = rmse( measured, func( xdata, solutionVectorChunks(:,pixel(1),pixel(2))));
    
    fprintf('\n');
    fprintf('Estimated T1 of pixel (%d/%d)\n', x, y);
    fprintf('  pixel-by-pixel, parfor:                    % 4.1f ms\n', solutionPixel(2, pixel(1), pixel(2)));
    fprintf('  vector, simultaneous fit:                  % 4.1f ms\n', solutionVector(2, pixel(1), pixel(2)));
    fprintf('  vector chunks, piecewise simultaneous fit: % 4.1f ms\n', solutionVectorChunks(2, pixel(1), pixel(2)));
    fprintf('\n');
    fprintf('RMSEs of pixel (%d/%d)\n', x, y);
    fprintf('  pixel-by-pixel, parfor:                    % 7.3f\n', rmsePixel);
    fprintf('  vector, simultaneous fit:                  % 7.3f\n', rmseVector);
    fprintf('  vector chunks, piecewise simultaneous fit: % 7.3f\n', rmseVectorChunks);

end
